#!/usr/bin/env python

from __future__ import division

import numpy as np
import os

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data

from wholecell.utils import units
from wholecell.utils.fitting import normalize, countsFromMassAndExpression

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

	# Load KB parameters

	bulkContainer = createBulkContainer(kb)

	# Mass
	
	g = growth_data.GrowthData(kb)
	massFractions60 = g.massFractions(60)
	rnaMass = massFractions60["rnaMass"]

	setRibosomeCountsConstrainedByPhysiology(kb, bulkContainer)
	

	rnaView = bulkContainer.countsView(kb.rnaData["id"])
	rnapView = bulkContainer.countsView(kb.rnapIds)
	proteinView = bulkContainer.countsView(kb.monomerData["id"])

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

	## Normalize expression and write out changes

	fitExpression(kb, bulkContainer)

	# Modify other properties

	fitRNAPolyTransitionRates(kb)

	## Calculate and set maintenance values

	# ----- Non growth associated maintenance -----
	kb.NGAM = NON_GROWTH_ASSOCIATED_MAINTENANCE * units.mmol / units.g / units.h

	# ----- Growth associated maintenance -----

	fitMaintenanceCosts(kb, bulkContainer)

# Sub-fitting functions

def createBulkContainer(kb):

	# Load from KB

	## IDs

	ids_molecules = kb.bulkMolecules['moleculeId']
	ids_rRNA23S = kb.rnaData["id"][kb.rnaData["isRRna23S"]]
	ids_rRNA16S = kb.rnaData["id"][kb.rnaData["isRRna16S"]]
	ids_rRNA5S = kb.rnaData["id"][kb.rnaData["isRRna5S"]]
	ids_tRNA = kb.rnaData["id"][kb.rnaData["isTRna"]]
	ids_mRNA = kb.rnaData["id"][kb.rnaData["isMRna"]]
	ids_protein = kb.monomerData["id"]

	## Mass fractions
	
	g = growth_data.GrowthData(kb)
	massFractions60 = g.massFractions(60)
	totalMass_RNA = massFractions60["rnaMass"]
	totalMass_protein = massFractions60["proteinMass"]

	totalMass_rRNA23S = totalMass_RNA * RRNA23S_MASS_SUB_FRACTION
	totalMass_rRNA16S = totalMass_RNA * RRNA16S_MASS_SUB_FRACTION
	totalMass_rRNA5S = totalMass_RNA * RRNA5S_MASS_SUB_FRACTION
	totalMass_tRNA = totalMass_RNA * TRNA_MASS_SUB_FRACTION
	totalMass_mRNA = totalMass_RNA * MRNA_MASS_SUB_FRACTION

	## Molecular weights

	individualMasses_rRNA23S = kb.getMass(ids_rRNA23S) / kb.nAvogadro
	individualMasses_rRNA16S = kb.getMass(ids_rRNA16S) / kb.nAvogadro
	individualMasses_rRNA5S = kb.getMass(ids_rRNA5S) / kb.nAvogadro
	individualMasses_tRNA = kb.rnaData["mw"][kb.rnaData["isTRna"]] / kb.nAvogadro
	individualMasses_mRNA = kb.rnaData["mw"][kb.rnaData["isMRna"]] / kb.nAvogadro
	individualMasses_protein = kb.monomerData["mw"] / kb.nAvogadro

	## Molecule distributions

	distribution_rRNA23S = np.array([1.] + [0.] * (ids_rRNA23S.size-1)) # currently only expressing first rRNA operon
	distribution_rRNA16S = np.array([1.] + [0.] * (ids_rRNA16S.size-1)) # currently only expressing first rRNA operon
	distribution_rRNA5S = np.array([1.] + [0.] * (ids_rRNA5S.size-1)) # currently only expressing first rRNA operon
	distribution_tRNA = normalize(kb.getTrnaAbundanceData(1 / units.h)['molar_ratio_to_16SrRNA'])
	distribution_mRNA = normalize(kb.rnaExpression['expression'][kb.rnaExpression['isMRna']])
	distribution_transcriptsByProtein = normalize(kb.rnaExpression['expression'][kb.rnaIndexToMonomerMapping])

	## Rates/times

	degradationRates = kb.monomerData["degRate"]
	doublingTime = kb.cellCycleLen

	# Construct bulk container

	bulkContainer = BulkObjectsContainer(ids_molecules, dtype = np.float64)

	## Assign rRNA counts based on mass

	totalCount_rRNA23S = totalCountFromMassesAndRatios(
		totalMass_rRNA23S,
		individualMasses_rRNA23S,
		distribution_rRNA23S
		)

	totalCount_rRNA23S.normalize()
	totalCount_rRNA23S.checkNoUnit()

	totalCount_rRNA16S = totalCountFromMassesAndRatios(
		totalMass_rRNA16S,
		individualMasses_rRNA16S,
		distribution_rRNA16S
		)

	totalCount_rRNA16S.normalize()
	totalCount_rRNA16S.checkNoUnit()

	totalCount_rRNA5S = totalCountFromMassesAndRatios(
		totalMass_rRNA5S,
		individualMasses_rRNA5S,
		distribution_rRNA5S
		)

	totalCount_rRNA5S.normalize()
	totalCount_rRNA5S.checkNoUnit()

	totalCount_rRNA_average = np.mean((totalCount_rRNA23S, totalCount_rRNA16S, totalCount_rRNA5S))

	counts_rRNA23S = totalCount_rRNA_average * distribution_rRNA23S
	counts_rRNA16S = totalCount_rRNA_average * distribution_rRNA16S
	counts_rRNA5S = totalCount_rRNA_average * distribution_rRNA5S

	bulkContainer.countsIs(counts_rRNA23S, ids_rRNA23S)
	bulkContainer.countsIs(counts_rRNA16S, ids_rRNA16S)
	bulkContainer.countsIs(counts_rRNA5S, ids_rRNA5S)

	## Assign tRNA counts based on mass and relative abundances (see Dong 1996)

	totalCount_tRNA = totalCountFromMassesAndRatios(
		totalMass_tRNA,
		individualMasses_tRNA,
		distribution_tRNA
		)

	totalCount_tRNA.normalize()
	totalCount_tRNA.checkNoUnit()

	counts_tRNA = totalCount_tRNA * distribution_tRNA

	bulkContainer.countsIs(counts_tRNA, ids_tRNA)

	## Assign mRNA counts based on mass and relative abundances (microarrays)

	totalCount_mRNA = totalCountFromMassesAndRatios(
		totalMass_mRNA,
		individualMasses_mRNA,
		distribution_mRNA
		)

	totalCount_mRNA.normalize()
	totalCount_mRNA.checkNoUnit()

	counts_mRNA = totalCount_mRNA * distribution_mRNA

	bulkContainer.countsIs(counts_mRNA, ids_mRNA)

	## Assign protein counts based on mass and mRNA counts
	distribution_protein = proteinDistribution(
		distribution_transcriptsByProtein,
		doublingTime,
		degradationRates
		)

	distribution_protein.checkNoUnit()

	distribution_protein = distribution_protein.asNumber() # remove units

	totalCount_protein = totalCountFromMassesAndRatios(
		totalMass_protein,
		individualMasses_protein,
		distribution_protein
		)

	totalCount_protein.normalize()
	totalCount_protein.checkNoUnit()
	
	counts_protein = totalCount_protein * distribution_protein

	bulkContainer.countsIs(counts_protein, ids_protein)

	return bulkContainer


def setRibosomeCountsConstrainedByPhysiology(kb, bulkContainer):
	'''
	setRibosomeCountsConstrainedByPhysiology

	Methodology: Set counts of ribosomal subunits based on three constraints.
	(1) Expected protein distribution doubles in one cell cycle
	(2) Measured rRNA mass fractions
	(3) Expected ribosomal subunit counts based on expression
	'''
	ribosome30SSubunits = kb.getComplexMonomers(kb.s30_fullComplex)[0]
	ribosome50SSubunits = kb.getComplexMonomers(kb.s50_fullComplex)[0]
	ribosome30SStoich = -1 * kb.getComplexMonomers(kb.s30_fullComplex)[1]
	ribosome50SStoich = -1 * kb.getComplexMonomers(kb.s50_fullComplex)[1]

	# -- CONSTRAINT 1: Expected protien distribution doubling -- #
	## Calculate minimium number of 30S and 50S subunits required in order to double our expected
	## protein distribution in one cell cycle
	proteinLengths = units.sum(kb.monomerData['aaCounts'], axis = 1)
	proteinDegradationRates =  kb.monomerData["degRate"]
	proteinCounts =  bulkContainer.counts(kb.monomerData["id"])

	nRibosomesNeeded = calculateMinPolymerizingEnzymeByProductDistribution(
	proteinLengths, kb.ribosomeElongationRate, proteinDegradationRates, proteinCounts, kb.cellCycleLen)
	nRibosomesNeeded.checkNoUnit()

	# Minimum number of ribosomes needed
	constraint1_ribosome30SCounts = (
		nRibosomesNeeded * ribosome30SStoich
		)

	constraint1_ribosome50SCounts = (
		nRibosomesNeeded * ribosome50SStoich
		)

	# -- CONSTRAINT 2: Measured rRNA mass fraction -- #
	## Calculate minimum number of 30S and 50S subunits based on measured mass fractions of
	## 16S, 23S, and 5S rRNA.
	rRna23SCounts = bulkContainer.counts(kb.rnaData["id"][kb.rnaData["isRRna23S"]])
	rRna16SCounts = bulkContainer.counts(kb.rnaData["id"][kb.rnaData["isRRna16S"]])
	rRna5SCounts = bulkContainer.counts(kb.rnaData["id"][kb.rnaData["isRRna5S"]])

	## 16S rRNA is in the 30S subunit
	massFracPredicted_30SCount = rRna16SCounts.sum()
	## 23S and 5S rRNA are in the 50S subunit
	massFracPredicted_50SCount = min(rRna23SCounts.sum(), rRna5SCounts.sum())

	constraint2_ribosome30SCounts = massFracPredicted_30SCount * ribosome30SStoich
	constraint2_ribosome50SCounts = massFracPredicted_50SCount * ribosome50SStoich

	# -- CONSTRAINT 3: Expected ribosomal subunit counts based distribution
	## Calculate fundamental ribosomal subunit count distribution based on RNA expression data
	## Already calculated and stored in bulkContainer
	ribosome30SCounts = bulkContainer.counts(ribosome30SSubunits)
	ribosome50SCounts = bulkContainer.counts(ribosome50SSubunits)

	# -- SET RIBOSOME FUNDAMENTAL SUBUNIT COUNTS TO MAXIMUM CONSTRAINT -- #
	bulkContainer.countsIs(
		np.fmax(np.fmax(ribosome30SCounts, constraint1_ribosome30SCounts), constraint2_ribosome30SCounts),
		ribosome30SSubunits
		)

	bulkContainer.countsIs(
		np.fmax(np.fmax(ribosome50SCounts, constraint1_ribosome50SCounts), constraint2_ribosome50SCounts),
		ribosome50SSubunits
		)
	
	# if np.any(ribosome30SView.counts() / ribosome30SStoich < nRibosomesNeeded) or np.any(ribosome50SView.counts() / ribosome50SStoich < nRibosomesNeeded):
	# 	raise NotImplementedError, "Cannot handle having too few ribosomes"


def calculateMinPolymerizingEnzymeByProductDistribution(productLengths, elongationRate, degradationRates, productCounts, cellCycleLen):
	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRate * (
			np.log(2) / cellCycleLen + degradationRates
			) * productCounts
		)
	return nPolymerizingEnzymeNeeded

def fitExpression(kb, bulkContainer):

	view_RNA = bulkContainer.countsView(kb.rnaData["id"])
	counts_protein = bulkContainer.counts(kb.monomerData["id"])

	g = growth_data.GrowthData(kb)
	massFractions60 = g.massFractions(60)
	totalMass_RNA = massFractions60["rnaMass"]

	### Modify kbFit to reflect our bulk container ###

	## RNA and monomer expression ##
	rnaExpressionContainer = BulkObjectsContainer(list(kb.rnaData["id"]), dtype = np.dtype("float64"))

	rnaExpressionContainer.countsIs(
		normalize(view_RNA.counts())
		)

	# Update mRNA expression to reflect monomer counts
	assert np.all(
		kb.monomerData["rnaId"][kb.monomerIndexToRnaMapping] == kb.rnaData["id"][kb.rnaData["isMRna"]]
		), "Cannot properly map monomer ids to RNA ids" # TODO: move to KB tests

	mRnaExpressionView = rnaExpressionContainer.countsView(kb.rnaData["id"][kb.rnaData["isMRna"]])
	mRnaExpressionFrac = np.sum(mRnaExpressionView.counts())

	mRnaExpressionView.countsIs(
		mRnaExpressionFrac * normalize(
			counts_protein *
			(np.log(2) / kb.cellCycleLen.asNumber(units.s) + kb.monomerData["degRate"].asNumber(1 / units.s))
			)[kb.monomerIndexToRnaMapping]
		)

	kb.rnaExpression['expression'] = rnaExpressionContainer.counts()

	# Set number of RNAs based on expression we just set
	nRnas = countsFromMassAndExpression(
		totalMass_RNA.asNumber(units.g),
		kb.rnaData["mw"].asNumber(units.g / units.mol),
		kb.rnaExpression['expression'],
		kb.nAvogadro.asNumber(1 / units.mol)
		)

	view_RNA.countsIs(nRnas * kb.rnaExpression['expression'])

	## Synthesis probabilities ##
	synthProb = normalize(
			(
			units.s * (
				np.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
				) * view_RNA.counts()
			).asNumber()
		)

	kb.rnaData["synthProb"][:] = synthProb


def fitRNAPolyTransitionRates(kb):
	## Transcription activation rate

	synthProb = kb.rnaData["synthProb"]
	rnaLengths = kb.rnaData["length"]

	elngRate = kb.rnaPolymeraseElongationRate

	# In our simplified model of RNA polymerase state transition, RNAp can be
	# active (transcribing) or inactive (free-floating).  To solve for the
	# rate of activation, we need to calculate the average rate of termination,
	# which is a function of the average transcript length and the 
	# transcription rate.

	averageTranscriptLength = units.dot(synthProb, rnaLengths)

	expectedTerminationRate = elngRate / averageTranscriptLength

	kb.transcriptionActivationRate = expectedTerminationRate * FRACTION_ACTIVE_RNAP / (1 - FRACTION_ACTIVE_RNAP)

	kb.fracActiveRnap = FRACTION_ACTIVE_RNAP


def fitMaintenanceCosts(kb, bulkContainer):
	aaCounts = kb.monomerData["aaCounts"]
	proteinCounts = bulkContainer.counts(kb.monomerData["id"])
	nAvogadro = kb.nAvogadro
	avgCellDryMassInit = kb.avgCellDryMassInit
	gtpPerTranslation = kb.gtpPerTranslation

	# GTPs used for translation (recycled, not incorporated into biomass)
	aaMmolPerGDCW = (
			units.sum(
				aaCounts * np.tile(proteinCounts.reshape(-1, 1), (1, 21)),
				axis = 0
			) * (
				(1 / (units.aa * nAvogadro)) *
				(1 / avgCellDryMassInit)
			)
		)

	aasUsedOverCellCycle = aaMmolPerGDCW.asNumber(units.mmol/units.g).sum()
	gtpUsedOverCellCycleMmolPerGDCW = gtpPerTranslation * aasUsedOverCellCycle

	darkATP = ( # This has everything we can't account for
		GROWTH_ASSOCIATED_MAINTENANCE -
		gtpUsedOverCellCycleMmolPerGDCW
		)

	# Assign the growth associated "dark energy" to translation
	# TODO: Distribute it amongst growth-related processes
	kb.gtpPerTranslation += darkATP / aasUsedOverCellCycle

# Math functions

def totalCountFromMassesAndRatios(totalMass, individualMasses, distribution):
	assert np.allclose(np.sum(distribution), 1)
	return 1 / units.dot(individualMasses, distribution) * totalMass


def proteinDistribution(distribution_mRNA, doublingTime, degradationRates):
	distributionUnnormed = 1 / (np.log(2) / doublingTime + degradationRates) * distribution_mRNA

	return distributionUnnormed / units.sum(distributionUnnormed)