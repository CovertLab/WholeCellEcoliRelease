#!/usr/bin/env python

from __future__ import division

import numpy as np
import os

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data

from wholecell.utils import units
from wholecell.utils.fitting import normalize

# Constants (should be moved to KB)
RRNA23S_MASS_SUB_FRACTION = 0.525 # This is the fraction of RNA that is 23S rRNA
RRNA16S_MASS_SUB_FRACTION = 0.271 # This is the fraction of RNA that is 16S rRNA
RRNA5S_MASS_SUB_FRACTION = 0.017 # This is the fraction of RNA that is 5S rRNA
TRNA_MASS_SUB_FRACTION = 0.146 # This is the fraction of RNA that is tRNA
MRNA_MASS_SUB_FRACTION = 0.041 # This is the fraction of RNA that is mRNA
GROWTH_ASSOCIATED_MAINTENANCE = 59.81 # mmol/gDCW (from Feist)
NON_GROWTH_ASSOCIATED_MAINTENANCE = 8.39 # mmol/gDCW/hr (from Feist)
FRACTION_ACTIVE_RNAP = 0.20 # from Dennis&Bremer; figure ranges from almost 100% to 20% depending on the growth rate

# Hacks
RNA_POLY_MRNA_DEG_RATE_PER_S = np.log(2) / 30. # half-life of 30 seconds
FRACTION_INCREASE_RIBOSOMAL_PROTEINS = 0.2  # reduce stochasticity from protein expression

# TODO: establish a controlled language for function behaviors (i.e. create* set* fit*)

FITNESS_THRESHOLD = 1e-9
MAX_FITTING_ITERATIONS = 100

def fitKb_1(kb):

	# Increase RNA poly mRNA deg rates
	# TODO: move to function
	# TODO: set this based on transcription unit structure
	# i.e. same synthesis prob. but different deg rates

	rnaPolySubunits = kb.process.complexation.getMonomers("APORNAP-CPLX[c]")["subunitIds"]

	subunitIndexes = np.array([np.where(kb.process.translation.monomerData["id"] == id_)[0].item() for id_ in rnaPolySubunits]) # there has to be a better way...

	mRNA_indexes = kb.relation.rnaIndexToMonomerMapping[subunitIndexes]

	kb.process.transcription.rnaData.struct_array["degRate"][mRNA_indexes] = RNA_POLY_MRNA_DEG_RATE_PER_S

	# nRNAs = kb.process.transcription.rnaData["expression"].size
	# kb.process.transcription.rnaData["expression"] = np.ones(nRNAs) / nRNAs
	# # WARNING - this doesn't update synthesis probabilities; for the moment, it doesn't need to

	# e = kb.process.transcription.rnaData["expression"].copy()

	# Fit synthesis probabilities for RNA

	for iteration in xrange(MAX_FITTING_ITERATIONS):

		initialExpression = kb.process.transcription.rnaData["expression"].copy()

		bulkContainer = createBulkContainer(kb)

		setRibosomeCountsConstrainedByPhysiology(kb, bulkContainer)

		setRNAPCountsConstrainedByPhysiology(kb, bulkContainer)

		# Normalize expression and write out changes

		fitExpression(kb, bulkContainer)

		finalExpression = kb.process.transcription.rnaData["expression"]

		degreeOfFit = np.sqrt(np.mean(np.square(initialExpression - finalExpression)))

		if degreeOfFit < FITNESS_THRESHOLD:
			break

	else:
		raise Exception("Fitting did not converge")

	# print np.linalg.norm(
	# 	e[kb.process.transcription.rnaData["isMRna"]] - kb.process.transcription.rnaData["expression"][kb.process.transcription.rnaData["isMRna"]],
	# 	2)

	# Modify other properties

	fitRNAPolyTransitionRates(kb, bulkContainer)

	## Calculate and set maintenance values

	# ----- Non growth associated maintenance -----
	kb.NGAM = NON_GROWTH_ASSOCIATED_MAINTENANCE * units.mmol / units.g / units.h

	# ----- Growth associated maintenance -----

	fitMaintenanceCosts(kb, bulkContainer)

# Sub-fitting functions

def createBulkContainer(kb):

	# Load from KB

	## IDs

	ids_molecules = kb.state.bulkMolecules.bulkData['id']
	ids_rRNA23S = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna23S"]]
	ids_rRNA16S = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna16S"]]
	ids_rRNA5S = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna5S"]]
	ids_tRNA = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isTRna"]]
	ids_mRNA = kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isMRna"]]
	ids_protein = kb.process.translation.monomerData["id"]

	## Mass fractions
	massFractions60 = kb.massFractions.massFractions(60 * units.min)

	totalMass_RNA = massFractions60["rnaMass"]
	totalMass_protein = massFractions60["proteinMass"]

	totalMass_rRNA23S = totalMass_RNA * RRNA23S_MASS_SUB_FRACTION
	totalMass_rRNA16S = totalMass_RNA * RRNA16S_MASS_SUB_FRACTION
	totalMass_rRNA5S = totalMass_RNA * RRNA5S_MASS_SUB_FRACTION
	totalMass_tRNA = totalMass_RNA * TRNA_MASS_SUB_FRACTION
	totalMass_mRNA = totalMass_RNA * MRNA_MASS_SUB_FRACTION

	## Molecular weights

	individualMasses_rRNA23S = kb.getter.getMass(ids_rRNA23S) / kb.constants.nAvogadro
	individualMasses_rRNA16S = kb.getter.getMass(ids_rRNA16S) / kb.constants.nAvogadro
	individualMasses_rRNA5S = kb.getter.getMass(ids_rRNA5S) / kb.constants.nAvogadro
	individualMasses_tRNA = kb.process.transcription.rnaData["mw"][kb.process.transcription.rnaData["isTRna"]] / kb.constants.nAvogadro
	individualMasses_mRNA = kb.process.transcription.rnaData["mw"][kb.process.transcription.rnaData["isMRna"]] / kb.constants.nAvogadro
	individualMasses_protein = kb.process.translation.monomerData["mw"] / kb.constants.nAvogadro

	## Molecule distributions

	distribution_rRNA23S = np.array([1.] + [0.] * (ids_rRNA23S.size-1)) # currently only expressing first rRNA operon
	distribution_rRNA16S = np.array([1.] + [0.] * (ids_rRNA16S.size-1)) # currently only expressing first rRNA operon
	distribution_rRNA5S = np.array([1.] + [0.] * (ids_rRNA5S.size-1)) # currently only expressing first rRNA operon
	distribution_tRNA = normalize(kb.massFractions.getTrnaAbundanceAtGrowthRate(60. * units.min)['molar_ratio_to_16SrRNA'])
	distribution_mRNA = normalize(kb.process.transcription.rnaData['expression'][kb.process.transcription.rnaData['isMRna']])
	distribution_transcriptsByProtein = normalize(kb.process.transcription.rnaData['expression'][kb.relation.rnaIndexToMonomerMapping])

	## Rates/times

	degradationRates = kb.process.translation.monomerData["degRate"]
	doublingTime = kb.constants.cellCycleLen

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

	totalCount_rRNA_average = sum([totalCount_rRNA23S, totalCount_rRNA16S, totalCount_rRNA5S]) / 3

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
	netLossRate_protein = netLossRateFromDilutionAndDegradation(doublingTime, degradationRates)

	distribution_protein = proteinDistributionFrommRNA(
		distribution_transcriptsByProtein,
		netLossRate_protein
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

	## HACK CODE: this code changes manually the number of EndoRNasess
	# endoRnaseIds = ["EG10856-MONOMER[p]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]"]
	# countsEndoRnase = bulkContainer.counts(endoRnaseIds)
	# normalizeDist = countsEndoRnase / countsEndoRnase.sum()
	# totalCounts = 6000.
	# newTotalCounts = totalCounts * normalizeDist
	# bulkContainer.countsIs(newTotalCounts, endoRnaseIds)
	## HACK CODE





	return bulkContainer


def setRibosomeCountsConstrainedByPhysiology(kb, bulkContainer):
	'''
	setRibosomeCountsConstrainedByPhysiology

	Methodology: Set counts of ribosomal subunits based on three constraints.
	(1) Expected protein distribution doubles in one cell cycle
	(2) Measured rRNA mass fractions
	(3) Expected ribosomal subunit counts based on expression
	'''
	ribosome30SSubunits = kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitIds']
	ribosome50SSubunits = kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitIds']
	ribosome30SStoich = kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitStoich']
	ribosome50SStoich = kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitStoich']

	# -- CONSTRAINT 1: Expected protien distribution doubling -- #
	## Calculate minimium number of 30S and 50S subunits required in order to double our expected
	## protein distribution in one cell cycle
	proteinLengths = units.sum(kb.process.translation.monomerData['aaCounts'], axis = 1)
	proteinDegradationRates =  kb.process.translation.monomerData["degRate"]
	proteinCounts =  bulkContainer.counts(kb.process.translation.monomerData["id"])

	netLossRate_protein = netLossRateFromDilutionAndDegradation(
		kb.constants.cellCycleLen,
		proteinDegradationRates
		)

	nRibosomesNeeded = calculateMinPolymerizingEnzymeByProductDistribution(
	proteinLengths, kb.constants.ribosomeElongationRate, netLossRate_protein, proteinCounts)
	nRibosomesNeeded.checkNoUnit()

	# Minimum number of ribosomes needed
	constraint1_ribosome30SCounts = (
		nRibosomesNeeded * ribosome30SStoich
		) * (1 + FRACTION_INCREASE_RIBOSOMAL_PROTEINS)

	constraint1_ribosome50SCounts = (
		nRibosomesNeeded * ribosome50SStoich
		) * (1 + FRACTION_INCREASE_RIBOSOMAL_PROTEINS)


	# -- CONSTRAINT 2: Measured rRNA mass fraction -- #
	## Calculate exact number of 30S and 50S subunits based on measured mass fractions of
	## 16S, 23S, and 5S rRNA.
	rRna23SCounts = bulkContainer.counts(kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna23S"]])
	rRna16SCounts = bulkContainer.counts(kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna16S"]])
	rRna5SCounts = bulkContainer.counts(kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna5S"]])

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

	# Fix rRNA counts
	bulkContainer.countsIs(rRna23SCounts, kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna23S"]])
	bulkContainer.countsIs(rRna16SCounts, kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna16S"]])
	bulkContainer.countsIs(rRna5SCounts, kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isRRna5S"]])


def setRNAPCountsConstrainedByPhysiology(kb, bulkContainer):
	# -- CONSTRAINT 1: Expected RNA distribution doubling -- #
	rnaLengths = units.sum(kb.process.transcription.rnaData['countsACGU'], axis = 1)
	rnaCounts = bulkContainer.counts(kb.process.transcription.rnaData['id'])
	## Number of endoRNases
	endoRnaseIds = ["EG10856-MONOMER[p]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]"]
	proteinCounts = bulkContainer.counts(kb.process.translation.monomerData["id"])
	endoRnaseCounts = np.zeros(len(endoRnaseIds))
	proteinIds = kb.process.translation.monomerData["id"]
	count_endoRN = 0

	for p in range (len(proteinIds)):
		if proteinIds[p] in endoRnaseIds:
			endoRnaseCounts[count_endoRN] = proteinCounts[p]
			count_endoRN += 1
	# print endoRnaseCounts
	print sum(endoRnaseCounts)



	# rnaLossRate = netLossRateFromDilutionAndDegradation(kb.cellCycleLen, kb.process.transcription.rnaData["degRate"])
	# import ipdb; ipdb.set_trace()
	rnaLossRate = netLossRateFromDilutionAndDegradationRNA(
		kb.constants.cellCycleLen,
		kb.process.transcription.rnaData["degRate"],
		#kb.constants.KcatEndoRNaseFullRNA, 
		kb.constants.KcatEndoRNasesFullRNA, 
		rnaCounts, 
		#np.sum(endoRnaseCounts)
		endoRnaseCounts
		)
	
	# nActiveRnapNeeded = calculateMinPolymerizingEnzymeByProductDistribution(rnaLengths, kb.rnaPolymeraseElongationRate, rnaLossRate, rnaCounts)
	nActiveRnapNeeded = calculateMinPolymerizingEnzymeByProductDistributionRNA(
		rnaLengths, kb.constants.rnaPolymeraseElongationRate, rnaLossRate)


	nActiveRnapNeeded.checkNoUnit()
	nRnapsNeeded = nActiveRnapNeeded / FRACTION_ACTIVE_RNAP

	minRnapSubunitCounts = (
		nRnapsNeeded * np.array([2, 1, 1, 1]) # Subunit stoichiometry # TODO: obtain automatically
		)

	# -- CONSTRAINT 2: Expected RNAP subunit counts based on distribution -- #
	rnapCounts = bulkContainer.counts(kb.moleculeGroups.rnapIds)

	## -- SET RNAP COUNTS TO MAXIMIM CONSTRAINTS -- #
	bulkContainer.countsIs(np.fmax(rnapCounts, minRnapSubunitCounts), kb.moleculeGroups.rnapIds)


def fitExpression(kb, bulkContainer):

	view_RNA = bulkContainer.countsView(kb.process.transcription.rnaData["id"])
	counts_protein = bulkContainer.counts(kb.process.translation.monomerData["id"])

	massFractions60 = kb.massFractions.massFractions(60 * units.min)
	totalMass_RNA = massFractions60["rnaMass"]

	doublingTime = kb.constants.cellCycleLen
	degradationRates_protein = kb.process.translation.monomerData["degRate"]

	netLossRate_protein = netLossRateFromDilutionAndDegradation(doublingTime, degradationRates_protein)

	### Modify kbFit to reflect our bulk container ###

	## RNA and monomer expression ##
	rnaExpressionContainer = BulkObjectsContainer(list(kb.process.transcription.rnaData["id"]), dtype = np.dtype("float64"))

	rnaExpressionContainer.countsIs(
		normalize(view_RNA.counts())
		)

	# Update mRNA expression to reflect monomer counts
	assert np.all(
		kb.process.translation.monomerData["rnaId"][kb.relation.monomerIndexToRnaMapping] == kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isMRna"]]
		), "Cannot properly map monomer ids to RNA ids" # TODO: move to KB tests

	mRnaExpressionView = rnaExpressionContainer.countsView(kb.process.transcription.rnaData["id"][kb.process.transcription.rnaData["isMRna"]])
	mRnaExpressionFrac = np.sum(mRnaExpressionView.counts())

	mRnaExpressionView.countsIs(
		mRnaExpressionFrac * mRNADistributionFromProtein(
			normalize(counts_protein), netLossRate_protein
			)[kb.relation.monomerIndexToRnaMapping]
		)

	kb.process.transcription.rnaData['expression'] = rnaExpressionContainer.counts()

	# Set number of RNAs based on expression we just set
	nRnas = totalCountFromMassesAndRatios(
		totalMass_RNA,
		kb.process.transcription.rnaData["mw"] / kb.constants.nAvogadro,
		kb.process.transcription.rnaData['expression']
		)

	nRnas.normalize()
	nRnas.checkNoUnit()

	view_RNA.countsIs(nRnas * kb.process.transcription.rnaData['expression'])

	## Number of endoRNases
	endoRnaseIds = ["EG10856-MONOMER[p]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]"]
	proteinCounts = counts_protein  # bulkContainer.counts(kb.process.translation.monomerDataData["id"])
	endoRnaseCounts = np.zeros(len(endoRnaseIds))
	proteinIds = kb.process.translation.monomerData["id"]
	count_endoRN = 0

	for p in range (len(proteinIds)):
		if proteinIds[p] in endoRnaseIds:
			endoRnaseCounts[count_endoRN] = proteinCounts[p]
			count_endoRN += 1
	#print endoRnaseCounts
	print np.sum(endoRnaseCounts)
	# import ipdb; ipdb.set_trace()

	## Synthesis probabilities ##
	# netLossRate_RNA = netLossRateFromDilutionAndDegradation(
	# 		doublingTime,
	# 		kb.process.transcription.rnaData["degRate"]
	# 	)
	# synthProb = normalize(
	# 		(
	# 		units.s
	# 			* netLossRate_RNA
	# 			* view_RNA.counts()
	# 		).asNumber()
	# 	)





	netLossRate_RNA = netLossRateFromDilutionAndDegradationRNA(
		(doublingTime / units.s).asNumber(),
		(kb.process.transcription.rnaData["degRate"] * units.s).asNumber(),
		#(kb.constants.KcatEndoRNaseFullRNA * units.s).asNumber(),
		(kb.constants.KcatEndoRNasesFullRNA * units.s).asNumber(),
		view_RNA.counts(), 
		#np.sum(endoRnaseCounts)
		endoRnaseCounts
		)
	# import ipdb; ipdb.set_trace()
	synthProb = normalize(netLossRate_RNA)

	kb.process.transcription.rnaData["synthProb"][:] = synthProb


def fitRNAPolyTransitionRates(kb, bulkContainer):
	## Transcription activation rate

	synthProb = kb.process.transcription.rnaData["synthProb"]

	rnaLengths = kb.process.transcription.rnaData["length"]


	elngRate = kb.constants.rnaPolymeraseElongationRate

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
	aaCounts = kb.process.translation.monomerData["aaCounts"]
	proteinCounts = bulkContainer.counts(kb.process.translation.monomerData["id"])
	nAvogadro = kb.constants.nAvogadro
	avgCellDryMassInit = kb.constants.avgCellDryMassInit
	gtpPerTranslation = kb.constants.gtpPerTranslation

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
	kb.constants.gtpPerTranslation += darkATP / aasUsedOverCellCycle


# Math functions

def totalCountFromMassesAndRatios(totalMass, individualMasses, distribution):
	"""
	Total mass = dot(mass, count)

	Fraction of i:
	f = count / Total counts

	Substituting:
	Total mass = dot(mass, f * Total counts)
	Total mass = Total counts * dot(mass, f)

	Total counts = Total mass / dot(mass, f)
	"""
	assert np.allclose(np.sum(distribution), 1)
	return 1 / units.dot(individualMasses, distribution) * totalMass


def proteinDistributionFrommRNA(distribution_mRNA, netLossRate):
	"""
	dP_i / dt = k * M_i - P_i * Loss_i

	At steady state:
	P_i = k * M_i / Loss_i

	Fraction of mRNA for ith gene is defined as:
	f_i = M_i / M_total

	Substituting in:
	P_i = k * f_i * M_total / Loss_i

	Normalizing P_i by summing over all i cancels out k and M_total
	assuming constant translation rate.
	"""

	assert np.allclose(np.sum(distribution_mRNA), 1)
	distributionUnnormed = 1 / netLossRate * distribution_mRNA

	return distributionUnnormed / units.sum(distributionUnnormed)


def mRNADistributionFromProtein(distribution_protein, netLossRate):
	"""
	dP_i / dt = k * M_i - P_i * Loss_i

	At steady state:
	M_i = Loss_i * P_i / k

	Fraction of protein for ith gene is defined as:
	f_i = P_i / P_total

	Substituting in:
	M_i = Loss_i * f_i * P_total / k

	Normalizing M_i by summing over all i cancles out k and P_total
	assuming a constant translation rate.

	"""
	assert np.allclose(np.sum(distribution_protein), 1)
	distributionUnnormed = netLossRate * distribution_protein

	return distributionUnnormed / units.sum(distributionUnnormed)


def calculateMinPolymerizingEnzymeByProductDistribution(productLengths, elongationRate, netLossRate, productCounts):
	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRate
			* netLossRate
			* productCounts
		)
	return nPolymerizingEnzymeNeeded

def calculateMinPolymerizingEnzymeByProductDistributionRNA(productLengths, elongationRate, netLossRate):
	nPolymerizingEnzymeNeeded = units.sum(
		productLengths / elongationRate
			* netLossRate
		)
	return nPolymerizingEnzymeNeeded


def netLossRateFromDilutionAndDegradation(doublingTime, degradationRates):
	return np.log(2) / doublingTime + degradationRates


def netLossRateFromDilutionAndDegradationRNA(doublingTime, degradationRates, kcatEndoRNase, RNACounts, endoRNaseCounts):
	# netRate = RNACounts * ( (np.log(2) / doublingTime) + degradationRates )

	# RNA decay considering RNase specificity (RNA lifetimes measured) and RNA accessibility (RNA counts)
	# new RNA decay model: kb - r ( ln(2)/tau + kcatEndoRN * EndoRN * kd / (Sum_g kd * r) ) = 0 
	# netRate = RNACounts * ( (np.log(2) / doublingTime) + (kcatEndoRNase * endoRNaseCounts * degradationRates / np.sum(degradationRates * RNACounts)) )
	# import ipdb; ipdb.set_trace()
	TotalEndoRNases = endoRNaseCounts * kcatEndoRNase
	netRate = RNACounts * ( (np.log(2) / doublingTime) + (TotalEndoRNases.sum() * degradationRates / np.sum(degradationRates * RNACounts)) )
	
	# RNA decay considering only RNase specificity (RNA lifetimes measured)
	# new RNA decay model: kb - kcatEndoRN * EndoRN * kd / (Sum_g kd * r) - r * ln(2) / tau = 0 
	# netRate = (sum(kcatEndoRNase * endoRNaseCounts) * degradationRates / np.sum(degradationRates)) + (np.log(2) / doublingTime * RNACounts)

	return netRate
	
