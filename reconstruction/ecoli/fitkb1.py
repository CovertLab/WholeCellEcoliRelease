#!/usr/bin/env python

from __future__ import division

import numpy as np
import os

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data

from wholecell.utils import units
from wholecell.utils.fitting import normalize, countsFromMassAndExpression, calcProteinCounts, calcProteinTotalCounts, calcProteinDistribution

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

def fitKb_1(kb):

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


def totalCountFromMassesAndRatios(totalMass, individualMasses, distribution):
	assert np.allclose(np.sum(distribution), 1)
	return totalMass / units.dot(individualMasses, distribution)


def setRRNACounts(kb, rnaMass, rRna23SView, rRna16SView, rRna5SView):

	## 23S rRNA Mass Fractions ##

	# Assume all 23S rRNAs are expressed equally
	nAvogadro = kb.nAvogadro

	nTypes_rRNA23S = kb.rnaData["isRRna23S"].sum()

	totalMass_rRNA23S = rnaMass * RRNA23S_MASS_SUB_FRACTION
	individualMasses_rRNA23S = kb.rnaData["mw"][kb.rnaData["isRRna23S"]] / nAvogadro # TODO: remove mws from rna data
	distribution_rRNA23S = np.array([1.] + [0.] * (nTypes_rRNA23S-1)) # currently only expressing first rRNA operon

	totalCount_rRNA23S = totalCountFromMassesAndRatios(
		totalMass_rRNA23S,
		individualMasses_rRNA23S,
		distribution_rRNA23S
		)

	totalCount_rRNA23S.checkNoUnit()

	## 16S rRNA Mass Fractions ##

	# Assume all 16S rRNAs are expressed equally
	nTypes_rRNA16S = kb.rnaData["isRRna16S"].sum()

	totalMass_rRNA16S = rnaMass * RRNA16S_MASS_SUB_FRACTION
	individualMasses_rRNA16S = kb.rnaData["mw"][kb.rnaData["isRRna16S"]] / nAvogadro # TODO: remove mws from rna data
	distribution_rRNA16S = np.array([1.] + [0.] * (nTypes_rRNA16S-1)) # currently only expressing first rRNA operon

	totalCount_rRNA16S = totalCountFromMassesAndRatios(
		totalMass_rRNA16S,
		individualMasses_rRNA16S,
		distribution_rRNA16S
		)

	totalCount_rRNA16S.checkNoUnit()

	## 5S rRNA Mass Fractions ##

	# Assume all 5S rRNAs are expressed equally

	nTypes_rRNA5S = kb.rnaData["isRRna5S"].sum()

	totalMass_rRNA5S = rnaMass * RRNA5S_MASS_SUB_FRACTION
	individualMasses_rRNA5S = kb.rnaData["mw"][kb.rnaData["isRRna5S"]] / nAvogadro # TODO: remove mws from rna data
	distribution_rRNA5S = np.array([1.] + [0.] * (nTypes_rRNA5S-1)) # currently only expressing first rRNA operon

	totalCount_rRNA5S = totalCountFromMassesAndRatios(
		totalMass_rRNA5S,
		individualMasses_rRNA5S,
		distribution_rRNA5S
		)

	totalCount_rRNA5S.checkNoUnit()

	# ## Correct numbers of 23S, 16S, 5S rRNAs so that they are all equal
	# # TODO: Maybe don't need to do this at some point (i.e., when the model is more sophisticated)
	totalCount_rRNA_average = np.mean((totalCount_rRNA23S, totalCount_rRNA16S, totalCount_rRNA5S))

	# Set counts

	rRna23SView.countsIs(totalCount_rRNA_average * distribution_rRNA23S)
	rRna16SView.countsIs(totalCount_rRNA_average * distribution_rRNA16S)
	rRna5SView.countsIs(totalCount_rRNA_average * distribution_rRNA5S)


def setTRNACounts(kb, rnaMass, tRnaView):

	nAvogadro = kb.nAvogadro

	## tRNA Mass Fractions ##

	# tRNA expression set based on data from Dong 1996
	totalMass_tRNA = rnaMass * TRNA_MASS_SUB_FRACTION
	individualMasses_tRNA = kb.rnaData["mw"][kb.rnaData["isTRna"]] / nAvogadro
	distribution_tRNA = normalize(kb.getTrnaAbundanceData(1 / units.h)['molar_ratio_to_16SrRNA'])

	totalCount_tRNA = totalCountFromMassesAndRatios(
		totalMass_tRNA,
		individualMasses_tRNA,
		distribution_tRNA
		)

	totalCount_tRNA.checkNoUnit()

	tRnaView.countsIs(totalCount_tRNA * distribution_tRNA)


def setMRNACounts(kb, rnaMass, mRnaView):

	nAvogadro = kb.nAvogadro

	## mRNA Mass Fractions ##
	totalMass_mRNA = rnaMass * MRNA_MASS_SUB_FRACTION
	individualMasses_mRNA = kb.rnaData["mw"][kb.rnaData["isMRna"]] / nAvogadro
	distribution_mRNA = normalize(kb.rnaExpression['expression'][kb.rnaExpression['isMRna']])

	totalCount_mRNA = totalCountFromMassesAndRatios(
		totalMass_mRNA,
		individualMasses_mRNA,
		distribution_mRNA
		)

	totalCount_mRNA.checkNoUnit()

	mRnaView.countsIs(totalCount_mRNA * distribution_mRNA)


def setMonomerCounts(kb, monomerMass, monomersView):

	monomersView.countsIs(calcProteinCounts(kb, monomerMass))