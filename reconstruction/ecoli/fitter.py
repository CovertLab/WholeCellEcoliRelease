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
from wholecell.utils.fitting import normalize, countsFromMassAndExpression, calcProteinCounts, calcProteinDistribution, calcProteinTotalCounts

from reconstruction.ecoli.fitkb1 import fitKb_1

def fitAtLevel(fitLevel, kb, simOutDir):
	# TODO: Obviously make this more sophisticated
	if fitLevel == 1:
		fitKb_1(kb)

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
