#!/usr/bin/env python

"""
Fitter

Adjust simulation parameters

TODO:
- document the math
- compute and use activation rates for RNA poly, ribosomes
- fit metabolism enzyme expression
- replace fake metabolite pools with measured metabolite pools

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/11/2013
"""

from __future__ import division

import numpy as np
import os
import collections

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data

from wholecell.utils import units
import unum # Imported here to be used in getCountsFromMassAndExpression assertions

# Constants (should be moved to KB)
RRNA23S_MASS_SUB_FRACTION = 0.525 # This is the fraction of RNA that is 23S rRNA
RRNA16S_MASS_SUB_FRACTION = 0.271 # This is the fraction of RNA that is 16S rRNA
RRNA5S_MASS_SUB_FRACTION = 0.017 # This is the fraction of RNA that is 5S rRNA
TRNA_MASS_SUB_FRACTION = 0.146 # This is the fraction of RNA that is tRNA
MRNA_MASS_SUB_FRACTION = 0.041 # This is the fraction of RNA that is mRNA
GROWTH_ASSOCIATED_MAINTENANCE = 59.81 # mmol/gDCW (from Feist)
NON_GROWTH_ASSOCIATED_MAINTENANCE = 8.39 # mmol/gDCW/hr (from Feist)
FRACTION_ACTIVE_RNAP = 0.20 # from Dennis&Bremer; figure ranges from almost 100% to 20% depending on the growth rate

# TODO: move many of these functions into another module


def fitAtLevel(fitLevel, kb, simOutDir):
	# TODO: Obviously make this more sophisticated
	if fitLevel == 1:
		fitKb(kb)

	if fitLevel == 2:
		# pass
		fitKb2(kb, simOutDir)

import tables

from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

N_SEEDS = 20

def fitKb2(kb, simOutDir):

	growthData = growth_data.GrowthData(kb)
	massFractions60 = growthData.massFractions(60)
	proteinMass = massFractions60["proteinMass"].asUnit(units.g)
	rnaMass = massFractions60["rnaMass"].asUnit(units.g)

	# Construct bulk container
	
	# We want to know something about the distribution of the copy numbers of 
	# macromolecules in the cell.  While RNA and protein expression can be
	# approximated using well-described statistical distributions, we need
	# absolute copy numbers to form complexes.  To get a distribution, we must
	# instantiate many cells, form complexes, and finally compute the 
	# statistics we will use in the fitting operations.

	bulkContainer = BulkObjectsContainer(kb.bulkMolecules['moleculeId'])
	rnaView = bulkContainer.countsView(kb.rnaData["id"])
	proteinView = bulkContainer.countsView(kb.monomerData["id"])
	complexationMoleculesView = bulkContainer.countsView(kb.complexationMoleculeNames)
	allMoleculesIDs = list(
		set(kb.rnaData["id"]) | set(kb.monomerData["id"]) | set(kb.complexationMoleculeNames)
		)
	allMoleculesView = bulkContainer.countsView(allMoleculesIDs)

	allMoleculeCounts = np.empty((N_SEEDS, allMoleculesView.counts().size), np.int64)

	complexationStoichMatrix = kb.complexationStoichMatrix().astype(np.int64, order = "F")

	complexationPrebuiltMatrices = mccBuildMatrices(
		complexationStoichMatrix
		)

	rnaDistribution = kb.rnaExpression['expression']

	rnaTotalCounts = countsFromMassAndExpression(
		rnaMass.asNumber(units.g),
		kb.rnaData["mw"].asNumber(units.g / units.mol),
		rnaDistribution,
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	proteinDistribution = calcProteinDistribution(kb)

	proteinTotalCounts = calcProteinTotalCounts(kb, proteinMass, proteinDistribution)
	
	for seed in xrange(N_SEEDS):
		randomState = np.random.RandomState(seed)

		allMoleculesView.countsIs(0)

		rnaView.countsIs(randomState.multinomial(
			rnaTotalCounts,
			rnaDistribution
			))

		proteinView.countsIs(randomState.multinomial(
			proteinTotalCounts,
			proteinDistribution
			))

		complexationMoleculeCounts = complexationMoleculesView.counts()

		updatedCompMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
			complexationMoleculeCounts,
			seed,
			complexationStoichMatrix,
			*complexationPrebuiltMatrices
			)

		complexationMoleculesView.countsIs(updatedCompMoleculeCounts)

		allMoleculeCounts[seed, :] = allMoleculesView.counts()

	bulkAverageContainer = BulkObjectsContainer(kb.bulkMolecules['moleculeId'], np.float64)
	bulkDeviationContainer = BulkObjectsContainer(kb.bulkMolecules['moleculeId'], np.float64)

	bulkAverageContainer.countsIs(allMoleculeCounts.mean(0), allMoleculesIDs)
	bulkDeviationContainer.countsIs(allMoleculeCounts.std(0), allMoleculesIDs)

	# Free up memory
	# TODO: make this more functional; one function for returning average & distribution
	del allMoleculeCounts
	del bulkContainer
	
	# ----- tRNA synthetase turnover rates ------
	# Fit tRNA synthetase kcat values based on expected rates of translation
	# compute values at initial time point

	## Compute rate of AA incorperation
	proteinComposition = kb.monomerData["aaCounts"]
	initialDryMass = kb.avgCellDryMassInit

	proteinMassFraction = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].asNumber(units.min) == 60.0
		]["proteinMassFraction"]

	initialProteinMass = initialDryMass * proteinMassFraction

	initialProteinCounts = calcProteinCounts(kb, initialProteinMass)

	initialProteinTranslationRate = (
		(np.log(2) / kb.cellCycleLen + kb.monomerData["degRate"]) * initialProteinCounts
		).asUnit(1 / units.s)

	initialAAPolymerizationRate = units.dot(
		units.transpose(proteinComposition), initialProteinTranslationRate
		).asUnit(units.aa / units.s)

	## Compute expression of tRNA synthetases
	## Assuming independence in variance
	synthetase_counts_by_group = np.zeros(len(kb.aa_synthetase_groups), dtype = np.float64)
	synthetase_variance_by_group = np.zeros(len(kb.aa_synthetase_groups), dtype = np.float)
	for idx, synthetase_group in enumerate(kb.aa_synthetase_groups.itervalues()):
		group_count = 0.
		group_variance = 0.
		for synthetase in synthetase_group:
			counts = bulkAverageContainer.countsView([synthetase]).counts()
			variance = bulkDeviationContainer.countsView([synthetase]).counts()
			group_count += counts
			group_variance += variance
		synthetase_counts_by_group[idx] = group_count
		synthetase_variance_by_group[idx] = group_variance
	
	## Saved for plotting
	kb.synthetase_counts = synthetase_counts_by_group
	kb.synthetase_variance = synthetase_variance_by_group
	kb.initial_aa_polymerization_rate = initialAAPolymerizationRate
	kb.minimum_trna_synthetase_rates = initialAAPolymerizationRate / synthetase_counts_by_group

	# TODO: Reimplement this with better fit taking into account the variance in aa
	#		utilization.
	## Scaling synthetase counts by -2*variance so that rates will be high enough
	## to accomodate stochastic behavior in the model without translation stalling.
	# scaled_synthetase_counts = synthetase_counts_by_group - (2 * synthetase_variance_by_group)
	scaled_synthetase_counts = synthetase_counts_by_group
	assert all(scaled_synthetase_counts > 0)

	predicted_trna_synthetase_rates = initialAAPolymerizationRate / scaled_synthetase_counts
	kb.trna_synthetase_rates = 2 * predicted_trna_synthetase_rates

	# fitKb2_metabolism(kb, simOutDir, bulkAverageContainer, bulkDeviationContainer)


from wholecell.utils.modular_fba import FluxBalanceAnalysis
def fitKb2_metabolism(kb, simOutDir, bulkAverageContainer, bulkDeviationContainer):

	# Load the simulation output

	## Effective biomass reaction
	with tables.open_file(os.path.join(simOutDir, "ConcentrationChange.hdf")) as h5file:
		time = h5file.root.ConcentrationChange.col("time")
		timeStep = h5file.root.ConcentrationChange.col("timeStep")

		# NOTE: units are M/s
		concentrationChange = h5file.root.ConcentrationChange.col("concentrationChange")

		names = h5file.root.names
		biomassMoleculeIDs = np.array(names.moleculeIDs.read())

	## Find the most extreme concentration flux, after removing the first few time steps

	# TODO: intelligent filtering - use the correlation coefficient?

	concentrationChange = concentrationChange[3:, :] # NOTE: magic number

	concentrationChangeMostPositive = concentrationChange.max(0)
	concentrationChangeMostNegative = concentrationChange.min(0)

	effectiveBiomassReaction = concentrationChangeMostPositive.copy()

	negativeIsMostExtreme = (np.abs(concentrationChangeMostNegative)
		> concentrationChangeMostPositive)

	effectiveBiomassReaction[negativeIsMostExtreme] = concentrationChangeMostNegative[negativeIsMostExtreme]

	## Build the enzyme-fitting problem

	# TODO: write a class for setting up LP problems

	values = []
	rowIndexes = []
	colIndexes = []

	lowerValues = []
	lowerIndexes = []

	upperValues = []
	upperIndexes = []

	objValues = []
	objIndexes = []

	rowNames = []
	colNames = []

	### Set up reverse reactions

	dt = 1 * units.s

	reactionStoich = kb.metabolismReactionStoich.copy()
	reactionEnzymes = kb.metabolismReactionEnzymes.copy()
	reactionRates = kb.metabolismReactionRates(dt)

	for reactionID in kb.metabolismReversibleReactions:
		reverseReactionID = "{} (reverse)".format(reactionID)
		assert reverseReactionID not in reactionStoich.viewkeys()
		reactionStoich[reverseReactionID] = {
			moleculeID:-coeff
			for moleculeID, coeff in reactionStoich[reactionID].viewitems()
			}

		if reactionID in reactionEnzymes.viewkeys():
			reactionEnzymes[reverseReactionID] = reactionEnzymes[reactionID]

		if reactionID in reactionRates.viewkeys():
			reactionRates[reverseReactionID] = reactionRates[reactionID]

	### Set up metabolites and biochemical reactions

	for reactionID, stoich in reactionStoich.viewitems():
		assert reactionID not in colNames
		reactionIndex = len(colNames)
		colNames.append(reactionID)

		for moleculeID, coeff in stoich.viewitems():
			try:
				moleculeIndex = rowNames.index(moleculeID)

			except ValueError:
				moleculeIndex = len(rowNames)
				rowNames.append(moleculeID)
			
			rowIndexes.append(moleculeIndex)
			colIndexes.append(reactionIndex)
			values.append(coeff)

	### Set up exchange reactions

	initWaterMass = kb.avgCellWaterMassInit
	initDryMass = kb.avgCellDryMassInit

	initCellMass = initWaterMass + initDryMass

	initCellVolume = initCellMass / kb.cellDensity

	coefficient = initDryMass / initCellVolume * dt

	exchangeConstraints = kb.metabolismExchangeConstraints(
		kb.metabolismExternalExchangeMolecules,
		coefficient,
		units.mmol / units.L
		)

	for moleculeID, constraint in zip(kb.metabolismExternalExchangeMolecules, exchangeConstraints):
		exchangeID = "{} exchange".format(moleculeID)

		assert exchangeID not in colNames
		exchangeIndex = len(colNames)
		colNames.append(exchangeID)

		moleculeIndex = rowNames.index(moleculeID)

		rowIndexes.append(moleculeIndex)
		colIndexes.append(exchangeIndex)
		values.append(-1.0)

		lowerIndexes.append(exchangeIndex)
		lowerValues.append(-min(constraint, 1e6))

	### Set up biomass reaction

	biomassID = "biomass reaction"
	assert biomassID not in colNames
	biomassIndex = len(colNames)
	colNames.append(biomassID)

	effectiveBiomassReaction *= 10**3 # change to mmol

	for moleculeID, coeff in zip(biomassMoleculeIDs, effectiveBiomassReaction):
		moleculeIndex = rowNames.index(moleculeID)

		rowIndexes.append(moleculeIndex)
		colIndexes.append(biomassIndex)
		values.append(-coeff)

		lowerIndexes.append(biomassIndex)
		lowerValues.append(+1) # must be capable of producing 100% of the biomass in a step

	### Set up enzyme usage

	enzymeRatesAll = collections.defaultdict(set)

	for reactionID, enzymeID in reactionEnzymes.viewitems():
		reactionRate = reactionRates[reactionID]

		enzymeRatesAll[enzymeID].add(reactionRate)

	enzymeIDs = enzymeRatesAll.keys()
	perEnzymeRates = {
		enzymeID:max(enzymeRates)
		for enzymeID, enzymeRates in enzymeRatesAll.viewitems()
		}

	minimalEnzymeCounts = np.fmax(
		bulkAverageContainer.counts(enzymeIDs) - 2 * bulkDeviationContainer.counts(enzymeIDs),
		0
		)

	enzymeConc = (
		1 / kb.nAvogadro / initCellVolume * minimalEnzymeCounts
		).asNumber(units.mmol / units.L)

	fullEnzymeRates = {
		enzymeID:perEnzymeRates[enzymeID] * enzymeConc[index]
		for index, enzymeID in enumerate(enzymeIDs)
		}

	for enzymeID, rateConstraint in fullEnzymeRates.viewitems():
		assert enzymeID not in rowNames
		enzymeIndex = len(rowNames)
		rowNames.append(enzymeID)

		constraintID = "{} constraint".format(enzymeID)
		assert constraintID not in colNames
		constraintIndex = len(colNames)
		colNames.append(constraintID)

		excessID = "{} excess capacity".format(enzymeID)
		assert excessID not in colNames
		excessIndex = len(colNames)
		colNames.append(excessID)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(constraintIndex)
		values.append(+1.0)

		upperIndexes.append(constraintIndex)
		upperValues.append(rateConstraint)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(excessIndex)
		values.append(+1.0)

		objIndexes.append(excessIndex)
		objValues.append(+1.0) # TODO: weighting

	for reactionID, enzymeID in reactionEnzymes.viewitems():
		if reactionID not in reactionRates.viewkeys():
			raise Exception("This code was not designed to handle enzymatic constraints without annotated rates.")

		reactionIndex = colNames.index(reactionID)
		enzymeIndex = rowNames.index(enzymeID)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(reactionIndex)
		values.append(-1)

	import cvxopt
	import cvxopt.solvers

	nRows = max(rowIndexes) + 1
	nCols = max(colIndexes) + 1

	assert len(values) == len(rowIndexes) == len(colIndexes)

	A = cvxopt.spmatrix(values, rowIndexes, colIndexes)

	b = cvxopt.matrix(np.zeros(nRows, np.float64))

	assert len(objIndexes) == len(objValues)

	objectiveFunction = np.zeros(nCols, np.float64)
	objectiveFunction[objIndexes] = objValues

	f = cvxopt.matrix(objectiveFunction)

	G = cvxopt.matrix(np.concatenate(
		[np.identity(nCols, np.float64), -np.identity(nCols, np.float64)]
		))

	assert len(upperIndexes) == len(upperValues)

	upperBound = np.empty(nCols, np.float64)
	upperBound.fill(1e6)
	upperBound[upperIndexes] = upperValues

	assert len(lowerIndexes) == len(lowerValues)

	lowerBound = np.empty(nCols, np.float64)
	lowerBound.fill(0)
	lowerBound[lowerIndexes] = lowerValues

	h = cvxopt.matrix(np.concatenate([upperBound, -lowerBound], axis = 0))

	oldOptions = cvxopt.solvers.options.copy()

	cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

	solution = cvxopt.solvers.lp(f, G, h, A, b, solver = "glpk")

	cvxopt.solvers.options.update(oldOptions)

	


def fitKb(kb):

	# Construct bulk container

	bulkContainer = BulkObjectsContainer(kb.bulkMolecules['moleculeId'], dtype = np.float64)

	rnaView = bulkContainer.countsView(kb.rnaData["id"])
	mRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isMRna"]])
	miscRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isMiscRna"]])
	rRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna"]])
	tRnaView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isTRna"]])

	rRna23SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna23S"]])
	rRna16SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna16S"]])
	rRna5SView = bulkContainer.countsView(kb.rnaData["id"][kb.rnaData["isRRna5S"]])

	monomersView = bulkContainer.countsView(kb.monomerData["id"])

	g = growth_data.GrowthData(kb)
	massFractions60 = g.massFractions(60)

	### RNA Mass fraction ###
	rnaMass = massFractions60["rnaMass"].asUnit(units.g)
	setRRNACounts(kb, rnaMass, rRna23SView, rRna16SView, rRna5SView)
	setTRNACounts(kb, rnaMass, tRnaView)
	setMRNACounts(kb, rnaMass, mRnaView)

	### Protein Mass fraction ###
	monomerMass = massFractions60["proteinMass"].asUnit(units.g)
	setMonomerCounts(kb, monomerMass, monomersView)

	### Ensure minimum numbers of enzymes critical for macromolecular synthesis ###

	rnapView = bulkContainer.countsView(kb.rnapIds)
	ribosome30SView = bulkContainer.countsView(kb.getComplexMonomers(kb.s30_fullComplex)[0])
	ribosome50SView = bulkContainer.countsView(kb.getComplexMonomers(kb.s50_fullComplex)[0])
	ribosome30SStoich = -1 * kb.getComplexMonomers(kb.s30_fullComplex)[1]
	ribosome50SStoich = -1 * kb.getComplexMonomers(kb.s50_fullComplex)[1]

	## Number of ribosomes needed ##
	monomerLengths = units.sum(kb.monomerData['aaCounts'], axis = 1)
	nRibosomesNeeded = units.sum(
		monomerLengths / kb.ribosomeElongationRate * (
			np.log(2) / kb.cellCycleLen + kb.monomerData["degRate"]
			) * monomersView.counts()
		).asNumber()
	
	# Minimum number of ribosomes needed
	sequencePredicted_min30SSubunitCounts = (
		nRibosomesNeeded * ribosome30SStoich
		)

	sequencePredicted_min50SSubunitCounts = (
		nRibosomesNeeded * ribosome50SStoich
		)

	# Number of ribosomes predicted from rRNA mass fractions
	# 16S rRNA is in the 30S subunit
	# 23S and 5S rRNA are in the 50S subunit
	massFracPredicted_30SCount = rRna16SView.counts().sum()
	massFracPredicted_50SCount = min(rRna23SView.counts().sum(), rRna5SView.counts().sum())
	massFracPrecicted_30SSubunitCounts = massFracPredicted_30SCount * ribosome30SStoich
	massFracPredicted_50SSubunitCounts = massFracPredicted_50SCount * ribosome50SStoich

	# Set ribosome subunit counts such that they are the maximum number from
	# (1) what is already in the container,
	# (2) what is predicted as needed based on sequence/elongation rate,
	# (3) what is predicted based on the rRNA mass fraction data
	ribosome30SView.countsIs(
		np.fmax(np.fmax(ribosome30SView.counts(), sequencePredicted_min30SSubunitCounts), massFracPrecicted_30SSubunitCounts)# + (1000 * ribosome30SStoich) # Added fudge factr of 1000
		)

	ribosome50SView.countsIs(
		np.fmax(np.fmax(ribosome50SView.counts(), sequencePredicted_min50SSubunitCounts), massFracPredicted_50SSubunitCounts)# + (1000 * ribosome50SStoich) # Added fudge factr of 1000
		)

	
	if np.any(ribosome30SView.counts() / ribosome30SStoich < nRibosomesNeeded) or np.any(ribosome50SView.counts() / ribosome50SStoich < nRibosomesNeeded):
		raise NotImplementedError, "Cannot handle having too few ribosomes"
	
	## Number of RNA Polymerases ##
	rnaLengths = units.sum(kb.rnaData['countsACGU'], axis = 1)

	nRnapsNeeded = units.sum(
		rnaLengths / kb.rnaPolymeraseElongationRate * (
			np.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
			) * rnaView.counts()
		).asNumber() / FRACTION_ACTIVE_RNAP

	minRnapCounts = (
		nRnapsNeeded * np.array([2, 1, 1, 1]) # Subunit stoichiometry # TODO: obtain automatically
		)

	rnapView.countsIs(
		np.fmax(rnapView.counts(), minRnapCounts)
		)

	### Modify kbFit to reflect our bulk container ###

	## RNA and monomer expression ##
	rnaExpressionContainer = BulkObjectsContainer(list(kb.rnaData["id"]), dtype = np.dtype("float64"))
	
	rnaExpressionContainer.countsIs(
		normalize(rnaView.counts())
		)

	# Update mRNA expression to reflect monomer counts
	assert np.all(
		kb.monomerData["rnaId"][kb.monomerIndexToRnaMapping] == kb.rnaData["id"][kb.rnaData["isMRna"]]
		), "Cannot properly map monomer ids to RNA ids"

	mRnaExpressionView = rnaExpressionContainer.countsView(kb.rnaData["id"][kb.rnaData["isMRna"]])
	mRnaExpressionFrac = np.sum(mRnaExpressionView.counts())

	mRnaExpressionView.countsIs(
		mRnaExpressionFrac * normalize(
			monomersView.counts() *
			(np.log(2) / kb.cellCycleLen.asNumber(units.s) + kb.monomerData["degRate"].asNumber(1 / units.s))
			)[kb.monomerIndexToRnaMapping]
		)

	kb.rnaExpression['expression'] = rnaExpressionContainer.counts()

	# Set number of RNAs based on expression we just set
	nRnas = countsFromMassAndExpression(
		rnaMass.asNumber(units.g),
		kb.rnaData["mw"].asNumber(units.g / units.mol),
		kb.rnaExpression['expression'],
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	rnaView.countsIs(nRnas * kb.rnaExpression['expression'])

	## Synthesis probabilities ##
	synthProb = normalize(
			(
			units.s * (
				np.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
				) * rnaView.counts()
			).asNumber()
		)

	kb.rnaData["synthProb"][:] = synthProb
	
	## Transcription activation rate

	# In our simplified model of RNA polymerase state transition, RNAp can be
	# active (transcribing) or inactive (free-floating).  To solve for the
	# rate of activation, we need to calculate the average rate of termination,
	# which is a function of the average transcript length and the 
	# transcription rate.

	averageTranscriptLength = units.dot(synthProb, rnaLengths)

	expectedTerminationRate = kb.rnaPolymeraseElongationRate / averageTranscriptLength

	kb.transcriptionActivationRate = expectedTerminationRate * FRACTION_ACTIVE_RNAP / (1 - FRACTION_ACTIVE_RNAP)

	kb.fracActiveRnap = FRACTION_ACTIVE_RNAP

	## Calculate and set maintenance values

	# ----- Non growth associated maintenance -----
	kb.NGAM = NON_GROWTH_ASSOCIATED_MAINTENANCE * units.mmol / units.g / units.h

	# ----- Growth associated maintenance -----

	# GTPs used for translation (recycled, not incorporated into biomass)
	aaMmolPerGDCW = (
			units.sum(
				kb.monomerData["aaCounts"] *
				np.tile(monomersView.counts().reshape(-1, 1), (1, 21)),
				axis = 0
			) * (
				(1 / (units.aa * kb.nAvogadro.asUnit(1 / units.mmol))) *
				(1 / kb.avgCellDryMassInit)
			)
		).asUnit(units.mmol / units.g)

	aasUsedOverCellCycle = aaMmolPerGDCW.asNumber().sum()
	gtpUsedOverCellCycleMmolPerGDCW = kb.gtpPerTranslation * aasUsedOverCellCycle

	darkATP = ( # This has everything we can't account for
		GROWTH_ASSOCIATED_MAINTENANCE -
		gtpUsedOverCellCycleMmolPerGDCW
		)

	# Assign the growth associated "dark energy" to translation
	# TODO: Distribute it amongst growth-related processes
	kb.gtpPerTranslation += darkATP / aaMmolPerGDCW.asNumber().sum()

def normalize(array):
	return np.array(array).astype("float") / np.linalg.norm(array, 1)

def countsFromMassAndExpression(mass, mws, relativeExpression, nAvogadro):
	"""
	countsFromMassAndExpression

	mass 				- float -				Total mass you want counts to sum to
	mws					- ndarray of floats -	Molecular weights of each species
	relativeExpression	- ndarray of floats	-	Relative expression of each species
	nAvogadro 			- float -				Avogadro's number

	Example:
		mass = 10.
		mws = [10., 5.]
		relativeExpression = [0.33, 0.66]
		countsFromMassAndExpression(mass, mws, relativeExpression, nAvogadro) = 1.93e23
	"""
	assert np.allclose(np.sum(relativeExpression), 1)
	assert type(mass) != unum.Unum
	assert type(mws) != unum.Unum
	assert type(relativeExpression) != unum.Unum
	assert type(nAvogadro) != unum.Unum
	return mass / np.dot(mws / nAvogadro, relativeExpression)

def setRRNACounts(kb, rnaMass, rRna23SView, rRna16SView, rRna5SView):

	## 23S rRNA Mass Fractions ##

	# Assume all 23S rRNAs are expressed equally
	rRna23SExpression = normalize(np.ones(rRna23SView.counts().size))

	nRRna23Ss = countsFromMassAndExpression(
		rnaMass.asNumber(units.g) * RRNA23S_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isRRna23S"]].asNumber(units.g / units.mol),
		rRna23SExpression,
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	## 16S rRNA Mass Fractions ##

	# Assume all 16S rRNAs are expressed equally
	rRna16SExpression = normalize(np.ones(rRna16SView.counts().size))

	nRRna16Ss = countsFromMassAndExpression(
		rnaMass.asNumber(units.g) * RRNA16S_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isRRna16S"]].asNumber(units.g / units.mol),
		rRna16SExpression,
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	## 5S rRNA Mass Fractions ##

	# Assume all 5S rRNAs are expressed equally
	rRna5SExpression = normalize(np.ones(rRna5SView.counts().size))

	nRRna5Ss = countsFromMassAndExpression(
		rnaMass.asNumber(units.g) * RRNA5S_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isRRna5S"]].asNumber(units.g / units.mol),
		rRna5SExpression,
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	# ## Correct numbers of 23S, 16S, 5S rRNAs so that they are all equal
	# # TODO: Maybe don't need to do this at some point (i.e., when the model is more sophisticated)
	nRRna23Ss = nRRna16Ss = nRRna5Ss = np.mean((nRRna23Ss, nRRna16Ss, nRRna5Ss))

	# TODO: Remove this hack once complexation is working
	rRna23SExpression[:] = 0.
	rRna23SExpression[0] = 1.

	rRna16SExpression[:] = 0.
	rRna16SExpression[0] = 1.

	rRna5SExpression[:] = 0.
	rRna5SExpression[0] = 1.

	rRna23SView.countsIs((nRRna23Ss * rRna23SExpression))
	rRna16SView.countsIs((nRRna16Ss * rRna16SExpression))
	rRna5SView.countsIs((nRRna5Ss * rRna5SExpression))


def setTRNACounts(kb, rnaMass, tRnaView):

	## tRNA Mass Fractions ##

	# tRNA expression set based on data from Dong 1996
	tRnaExpression = normalize(kb.getTrnaAbundanceData(1 / units.h)['molar_ratio_to_16SrRNA'])

	nTRnas = countsFromMassAndExpression(
		rnaMass.asNumber(units.g) * TRNA_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isTRna"]].asNumber(units.g / units.mol),
		tRnaExpression,
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	tRnaView.countsIs((nTRnas * tRnaExpression))

def setMRNACounts(kb, rnaMass, mRnaView):

	## mRNA Mass Fractions ##

	mRnaExpression = normalize(kb.rnaExpression['expression'][kb.rnaExpression['isMRna']])

	nMRnas = countsFromMassAndExpression(
		rnaMass.asNumber(units.g) * MRNA_MASS_SUB_FRACTION,
		kb.rnaData["mw"][kb.rnaData["isMRna"]].asNumber(units.g / units.mol),
		mRnaExpression,
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	mRnaView.countsIs((nMRnas * mRnaExpression))


def setMonomerCounts(kb, monomerMass, monomersView):

	monomersView.countsIs(calcProteinCounts(kb, monomerMass))


def calcProteinCounts(kb, monomerMass):
	monomerExpression = calcProteinDistribution(kb)

	nMonomers = calcProteinTotalCounts(kb, monomerMass, monomerExpression)

	return nMonomers * monomerExpression


def calcProteinTotalCounts(kb, monomerMass, monomerExpression):
	return countsFromMassAndExpression(
		monomerMass.asNumber(units.g),
		kb.monomerData["mw"].asNumber(units.g / units.mol),
		monomerExpression,
		kb.nAvogadro.asNumber(1 / units.mol)
		)

def calcProteinDistribution(kb):
	return normalize(
		kb.rnaExpression['expression'][kb.rnaIndexToMonomerMapping] /
		(np.log(2) / kb.cellCycleLen.asNumber(units.s) + kb.monomerData["degRate"].asNumber(1 / units.s))
		)


if __name__ == "__main__":
	import wholecell.utils.constants

	kb = cPickle.load(
		open(os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
			), "rb")
		)
	
	fitKb(kb)
