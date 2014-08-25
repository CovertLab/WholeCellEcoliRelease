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
import copy
import collections

import wholecell.states.bulk_molecules
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data

from wholecell.utils import units
import unum

# Constants (should be moved to KB)
RRNA23S_MASS_SUB_FRACTION = 0.525 # This is the fraction of RNA that is 23S rRNA
RRNA16S_MASS_SUB_FRACTION = 0.271 # This is the fraction of RNA that is 16S rRNA
RRNA5S_MASS_SUB_FRACTION = 0.017 # This is the fraction of RNA that is 5S rRNA
TRNA_MASS_SUB_FRACTION = 0.146 # This is the fraction of RNA that is tRNA
MRNA_MASS_SUB_FRACTION = 0.041 # This is the fraction of RNA that is mRNA
GROWTH_ASSOCIATED_MAINTENANCE = 59.81 # mmol/gDCW (from Feist)
NON_GROWTH_ASSOCIATED_MAINTENANCE = 8.39 # mmol/gDCW/hr (from Feist)

# Correction factors
EXCESS_RNAP_CAPACITY = 1.
EXCESS_FREE_DNTP_CAPACITY = 1.3
# If RNA-poly capacity exactly matches the amount needed to double RNAs over a 
# cell cycle, the simulation will be unable to double RNAs since a small number
# of RNA-polymerases must be turned over following termination.  It may be 
# possible to choose an excess capacity coefficient rationally based on 
# diffusive limitations, i.e., one that does not depend on simulation 
# particulars, but this has yet to be explored.

# Fitter logic 
# TODO: confirm this with Derek
# TODO: split off these subroutines in the main fitter function
# 1) Assign expected quantities based on dry mass composition, expression, and sequences
# 2) Ensure that there is enough RNAP/ribosome capacity for (1), and adjust if needed
# 3) Update the metabolism FBA objective based on expression

# TODO: move many of these functions into another module


def fitAtLevel(fitLevel, kb, simOutDir):
	# TODO: Obviously make this more sophisticated
	if fitLevel == 1:
		fitKb(kb)

	if fitLevel == 2:
		pass
		# fitKb2(kb, simOutDir)

import tables

from wholecell.utils.modular_fba import FluxBalanceAnalysis

def fitKb2(kb, simOutDir):

	# Load the simulation output

	## Effective biomass reactuib
	with tables.open_file(os.path.join(simOutDir, "ConcentrationChange.hdf")) as h5file:
		time = h5file.root.ConcentrationChange.col("time")
		timeStep = h5file.root.ConcentrationChange.col("timeStep")

		# NOTE: units are M/s
		concentrationChange = h5file.root.ConcentrationChange.col("concentrationChange")

		names = h5file.root.names
		moleculeIDs = np.array(names.moleculeIDs.read())

	## Find the most extreme concentration flux, after removing the first few time steps

	# TODO: intelligent filtering - use the correlation coefficient?

	concentrationChange = concentrationChange[3:, :] # NOTE: magic number

	concentrationChangeMostPositive = concentrationChange.max(0)
	concentrationChangeMostNegative = concentrationChange.min(0)

	effectiveBiomassReaction = concentrationChangeMostPositive.copy()

	negativeIsMostExtreme = (np.abs(concentrationChangeMostNegative)
		> concentrationChangeMostPositive)

	effectiveBiomassReaction[negativeIsMostExtreme] = concentrationChangeMostNegative[negativeIsMostExtreme]

	## Build the standard FBA problem and the MOMA problem

	DELTA_T = 1 # use a time-step of one-second (doesn't really matter)

	# objective = dict(zip(moleculeIDs, effectiveBiomassReaction*DELTA_T))
	objective = dict(zip(moleculeIDs, effectiveBiomassReaction*DELTA_T*10**3))

	reactionRates = kb.metabolismReactionRates(DELTA_T * units.s)
	for reactionID, reactionRate in reactionRates.viewitems():
		reactionRates[reactionID] = max(reactionRate, 500) # TODO: remove this hack

	fba = FluxBalanceAnalysis(
		kb.metabolismReactionStoich.copy(),
		kb.metabolismExternalExchangeMolecules,
		objective,
		objectiveType = "standard",
		reversibleReactions = kb.metabolismReversibleReactions,
		reactionEnzymes = kb.metabolismReactionEnzymes.copy(), # TODO: copy in class
		reactionRates = reactionRates,
		)

	externalMoleculeIDs = fba.externalMoleculeIDs()

	initWaterMass = kb.avgCellWaterMassInit
	initDryMass = kb.avgCellDryMassInit

	initCellMass = (
		initWaterMass
		+ initDryMass
		)

	initCellVolume = initCellMass / kb.cellDensity

	coefficient = initDryMass / initCellVolume * (DELTA_T * units.s)

	externalMoleculeLevels = kb.metabolismExchangeConstraints(
		externalMoleculeIDs,
		coefficient,
		# units.mol / units.L
		units.mmol / units.L
		)

	fba.externalMoleculeLevelsIs(externalMoleculeLevels)

	## Calculate protein concentrations and assign enzymatic limits
	growthData = growth_data.GrowthData(kb)
	massFractions60 = growthData.massFractions(60)
	proteinMass = massFractions60["proteinMass"].asUnit(units.g)

	proteinConc = dict(zip(
		kb.monomerData["id"],
		calcProteinCounts(kb, proteinMass) / kb.nAvogadro.asNumber(1 / units.mmol) / initCellVolume.asNumber(units.L)
		))

	enzymeConc = np.array([
		proteinConc[enzymeID] if enzymeID in proteinConc.viewkeys() else np.inf
		for enzymeID in fba.enzymeIDs()
		])

	unaccountedEnzymes = [
		enzymeID for enzymeID in fba.enzymeIDs() if enzymeID not in proteinConc
		] # these should all be complexes, need to handle

	# TODO: compute complex concentrations

	# TODO: scaling factor in FBA solver

	fba.enzymeLevelsIs(enzymeConc)

	fba.run()

	print fba.objectiveReactionFlux()

	import ipdb; ipdb.set_trace()


def fitKb(kb):

	# Construct bulk container

	bulkContainer = wholecell.states.bulk_molecules.bulkObjectsContainer(kb, dtype = np.dtype("float64"))

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
	setRNACounts(
		kb, rnaMass, mRnaView,
		rRna23SView, rRna16SView, rRna5SView, tRnaView
		)


	### Protein Mass fraction ###
	monomerMass = massFractions60["proteinMass"].asUnit(units.g)
	setMonomerCounts(kb, monomerMass, monomersView)

	### Ensure minimum numbers of enzymes critical for macromolecular synthesis ###

	rnapView = bulkContainer.countsView(kb.rnapIds)

	## Number of ribosomes needed ##
	monomerLengths = units.sum(kb.monomerData['aaCounts'], axis = 1)
	nRibosomesNeeded = units.sum(
		monomerLengths / kb.ribosomeElongationRate * (
			np.log(2) / kb.cellCycleLen + kb.monomerData["degRate"]
			) * monomersView.counts()
		).asNumber()
	
	if np.sum(rRna23SView.counts()) < nRibosomesNeeded:
		raise NotImplementedError, "Cannot handle having too few ribosomes"

	## Number of RNA Polymerases ##
	rnaLengths = units.sum(kb.rnaData['countsACGU'], axis = 1)

	nRnapsNeeded = units.sum(
		rnaLengths / kb.rnaPolymeraseElongationRate * (
			np.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
			) * rnaView.counts()
		).asNumber() * EXCESS_RNAP_CAPACITY

	minRnapCounts = (
		nRnapsNeeded * np.array([2, 1, 1, 1]) # Subunit stoichiometry
		)

	rnapView.countsIs(
		np.fmax(rnapView.counts(), minRnapCounts)
		)

	
	### Modify kbFit to reflect our bulk container ###

	## Fraction of active Ribosomes ##
	kb.fracActiveRibosomes = float(nRibosomesNeeded) / np.sum(rRna23SView.counts())

	## RNA and monomer expression ##
	rnaExpressionContainer = wholecell.containers.bulk_objects_container.BulkObjectsContainer(list(kb.rnaData["id"]), dtype = np.dtype("float64"))

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
		mRnaExpressionFrac * normalize(monomersView.counts()[kb.monomerIndexToRnaMapping])
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


	# ----- tRNA synthetase ------
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
	synthetase_counts_by_group = np.zeros(len(kb.aa_synthetase_groups), dtype = np.float64)
	for idx, synthetase_group in enumerate(kb.aa_synthetase_groups.itervalues()):
		group_count = 0.
		for synthetase in synthetase_group:
			try:
				subunits, stoich = kb.getComplexMonomers(synthetase)
				subunitCounts = bulkContainer.countsView(subunits).counts()
				group_count += np.min(subunitCounts / (-1*stoich))
			except IndexError:
				counts = bulkContainer.countsView([synthetase]).counts()
				group_count += counts
			synthetase_counts_by_group[idx] = group_count
	
	predicted_trna_synthetase_rates = initialAAPolymerizationRate / synthetase_counts_by_group
	kb.trna_synthetase_rates = predicted_trna_synthetase_rates

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

def setRNACounts(kb, rnaMass, mRnaView, rRna23SView, rRna16SView, rRna5SView, tRnaView):

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
	monomerExpression = normalize(
		kb.rnaExpression['expression'][kb.rnaIndexToMonomerMapping] /
		(np.log(2) / kb.cellCycleLen.asNumber(units.s) + kb.monomerData["degRate"].asNumber(1 / units.s))
		)

	nMonomers = countsFromMassAndExpression(
		monomerMass.asNumber(units.g),
		kb.monomerData["mw"].asNumber(units.g / units.mol),
		monomerExpression,
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	return nMonomers * monomerExpression


if __name__ == "__main__":
	import wholecell.utils.constants

	kb = cPickle.load(
		open(os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
			), "rb")
		)
	
	fitKb(kb)
