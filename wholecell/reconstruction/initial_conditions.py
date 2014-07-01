
from __future__ import division

from itertools import izip

import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.reconstruction.fitter import countsFromMassAndExpression
from wholecell.reconstruction.fitter import normalize

def calcInitialConditions(sim, kb):
	randomState = sim.randomState

	timeStep = sim.timeStepSec # This is a poor solution but will suffice for now

	bulkMolCntr = sim.states['BulkMolecules'].container
	uniqueMolCntr = sim.states["UniqueMolecules"].container
	bulkChrmCntr = sim.states["BulkChromosome"].container

	# Set up states
	initializeBulkMolecules(bulkMolCntr, kb, randomState, timeStep)
	initializeBulkChromosome(bulkChrmCntr, kb, randomState, timeStep)
	initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, kb, randomState, timeStep)


def initializeBulkMolecules(bulkMolCntr, kb, randomState, timeStep):

	## Set protein counts from expression
	initializeProteinMonomers(bulkMolCntr, kb, randomState, timeStep)

	## Set RNA counts from expression
	initializeRNA(bulkMolCntr, kb, randomState, timeStep)

	## Set DNA
	initializeDNA(bulkMolCntr, kb, randomState, timeStep)

	## Set other biomass components
	initializeBulkComponents(bulkMolCntr, kb, randomState, timeStep)

	## Set pools
	initializePools(bulkMolCntr, kb, randomState, timeStep)

	## Set water
	initializeBulkWater(bulkMolCntr, kb, randomState, timeStep)

	## Form complexes
	initializeComplexes(bulkMolCntr, kb, randomState, timeStep)


def initializeBulkChromosome(bulkChrmCntr, kb, randomState, timeStep):
	## Set genes
	initializeGenes(bulkChrmCntr, kb, timeStep)


def initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, kb, randomState, timeStep):
	initializeTranscription(bulkMolCntr, uniqueMolCntr, kb, randomState, timeStep)
	initializeTranslation(bulkMolCntr, uniqueMolCntr, kb, randomState, timeStep)
	initializeReplication(uniqueMolCntr, kb)


def initializeProteinMonomers(bulkMolCntr, kb, randomState, timeStep):
	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60
		]

	monomersView = bulkMolCntr.countsView(kb.monomerData["id"])
	monomerMassFraction = float(dryComposition60min["proteinMassFraction"])
	monomerMass = kb.avgCellDryMassInit.to('DCW_g') * monomerMassFraction

	monomerExpression = normalize(kb.rnaExpression['expression'][kb.rnaIndexToMonomerMapping].to('dimensionless'))

	nMonomers = countsFromMassAndExpression(
		monomerMass.to("DCW_g").magnitude,
		kb.monomerData["mw"].to("g/mol").magnitude,
		monomerExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	monomersView.countsIs(
		randomState.multinomial(nMonomers, monomerExpression)
		)

	# monomersView.countsIs(nMonomers * monomerExpression)


def initializeRNA(bulkMolCntr, kb, randomState, timeStep):
	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].magnitude == 60
		]

	rnaView = bulkMolCntr.countsView(kb.rnaData["id"])
	rnaMassFraction = float(dryComposition60min["rnaMassFraction"])
	rnaMass = kb.avgCellDryMassInit * rnaMassFraction

	rnaExpression = normalize(kb.rnaExpression['expression'].to('dimensionless'))

	nRnas = countsFromMassAndExpression(
		rnaMass.to("DCW_g").magnitude,
		kb.rnaData["mw"].to('g/mol').magnitude,
		rnaExpression,
		kb.nAvogadro.to('1/mol').magnitude
		)

	rnaView.countsIs(
		randomState.multinomial(nRnas, rnaExpression)
		)

	# rnaView.countsIs(nRnas * rnaExpression)


def initializeDNA(bulkMolCntr, kb, randomState, timeStep):

	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].magnitude == 60
		]

	dntpsView = bulkMolCntr.countsView(kb.dNtpIds)
	dnmpsView = bulkMolCntr.countsView(kb.dNmpNuclearIds)
	dnaMassFraction = float(dryComposition60min["dnaMassFraction"])
	dnaMass = kb.avgCellDryMassInit.magnitude * dnaMassFraction

	dnaExpression = normalize(np.array([
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		], dtype = np.float64))

	dnmpMws = kb.getMass(kb.dNmpNuclearIds).magnitude

	dntpMws = kb.getMass(kb.dNtpIds).magnitude

	dnmpsView.countsIs([
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		])

	chromMass = (
		np.dot(dnmpsView.counts(), dnmpMws) - 2 * kb.genomeLength * 17.01
		) / kb.nAvogadro.magnitude

	# dNTP Pools are handled here
	nDntps = countsFromMassAndExpression(
		dnaMass - chromMass,
		dntpMws,
		dnaExpression,
		kb.nAvogadro.magnitude
		)
	dntpsView.countsIs(
		randomState.multinomial(nDntps, dnaExpression)
		)


def initializeBulkComponents(bulkMolCntr, kb, randomState, timeStep):
	biomassContainer = BulkObjectsContainer(
		list(kb.wildtypeBiomass["metaboliteId"]), dtype = np.dtype("float64")
		)
	biomassContainer.countsIs(
		kb.wildtypeBiomass["biomassFlux"].to("millimole/DCW_gram").magnitude
		)

	notPRDMetabolites = (
		list(kb.cellGlycogenFractionData["metaboliteId"]) +
		list(kb.cellMureinFractionData["metaboliteId"]) +
		list(kb.cellLPSFractionData["metaboliteId"]) +
		list(kb.cellLipidFractionData["metaboliteId"]) +
		list(kb.cellInorganicIonFractionData["metaboliteId"]) +
		list(kb.cellSolublePoolFractionData["metaboliteId"])
		)

	notPRDBulkView = bulkMolCntr.countsView(notPRDMetabolites)

	notPRDBiomassView = biomassContainer.countsView(notPRDMetabolites)

	notPRDBulkView.countsIs((
		kb.avgCellDryMassInit.to("DCW_gram").magnitude *
		notPRDBiomassView.counts() *
		kb.nAvogadro.to("1 / millimole").magnitude
		))


def initializePools(bulkMolCntr, kb, randomState, timeStep):
	# Note: This is adding dry biomass, so the cell will appear heavier

	from wholecell.reconstruction.knowledge_base_ecoli import AMINO_ACID_1_TO_3_ORDERED

	biomassContainer = BulkObjectsContainer(
		list(kb.wildtypeBiomass["metaboliteId"]), dtype = np.dtype("float64")
		)
	biomassContainer.countsIs(
		kb.wildtypeBiomass["biomassFlux"].to("millimole/DCW_gram").magnitude
		)

	ntpsBiomassView = biomassContainer.countsView(kb.ntpIds)
	dntpsBiomassView = biomassContainer.countsView(kb.dNtpIds)
	aasBiomassView = biomassContainer.countsView(kb.aaIDs)

	ppiBulkView = bulkMolCntr.countView("PPI[c]")
	ntpsBulkView = bulkMolCntr.countsView(kb.ntpIds)
	dntpsBulkView = bulkMolCntr.countsView(kb.dNtpIds)
	aasBulkView = bulkMolCntr.countsView(kb.aaIDs)

	dt = timeStep
	tau_d = kb.cellCycleLen.to("second").magnitude

	##### dNTPs are handled in DNA initialization because they are special
	## dNTPs
	# dntpsFromLastStep = (
	# 	dntpsBiomassView.counts() * (1 - np.exp(-np.log(2) / tau_d * dt)) *
	# 	kb.nAvogadro.to("1 / millimole").magnitude *
	# 	kb.avgCellDryMassInit.to("DCW_gram").magnitude
	# 	)
	# dntpsBulkView.countsIs(dntpsFromLastStep)

	## NTPs (for Transcription)
	ntpsFromLastStep = (
		ntpsBiomassView.counts() *
		(1 - np.exp(-np.log(2) / tau_d * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		).astype(np.int64)
	ntpsBulkView.countsInc(ntpsFromLastStep)

	## Small molecules for Translation
	gtpInitial = kb.gtpPoolSize * dt * kb.nAvogadro.to("1/millimole") * kb.avgCellDryMassInit.to("DCW_gram")
	bulkMolCntr.countsInc(
		gtpInitial.magnitude.astype(np.int64),
		["GTP[c]", "PPI[c]", "GMP[c]"]
		)

	## Small molecules for GAM
	atpInitial = kb.atpPoolSize * dt * kb.nAvogadro.to("1/millimole") * kb.avgCellDryMassInit.to("DCW_gram")
	bulkMolCntr.countsInc(
		atpInitial.magnitude.astype(np.int64),
		["ATP[c]", "H2O[c]", "PI[c]", "ADP[c]", "H[c]"]
		)

	## Amino Acids
	aasFromLastStep = (
		aasBiomassView.counts() * (1 - np.exp(-np.log(2) / tau_d * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)
	aasBulkView.countsIs(aasFromLastStep)

	## PPI
	# Assumption:
	# NTPs and dNTPs from two steps ago got completely used up in the last step
	ntpsFromTwoStepsAgo = (
		ntpsBiomassView.counts() *
		(np.exp(-np.log(2) / tau_d * dt) - np.exp(-np.log(2) / tau_d * 2 * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)

	# TODO: This is definitely wrong because replication isn't metabolically limited
	dntpsFromTwoStepsAgo = (
		dntpsBiomassView.counts() *
		(np.exp(-np.log(2) / tau_d * dt) - np.exp(-np.log(2) / tau_d * 2 * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)

	ppiFromNtps = np.round(np.sum(
		ntpsFromLastStep
		))

	# ppiFromDntps = np.round(np.sum(
	# 	dntpsFromLastStep
	# 	))
	ppiFromDntps = 0

	ppiBulkView.countInc(ppiFromNtps + ppiFromDntps)


def initializeBulkWater(bulkMolCntr, kb, randomState, timeStep):
	h2oView = bulkMolCntr.countView('H2O[c]')

	nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
	mwH2O = kb.getMass(['H2O'])[0]
	avgCellWaterMassInit = kb.avgCellWaterMassInit.to('water_g').magnitude

	h2oView.countIs(
		(avgCellWaterMassInit) / mwH2O * nAvogadro
		)


def initializeGenes(bulkChrmCntr, kb, timeStep):
	"""
	initializeGenes

	Purpose:
	Initalizes the counts of genes in BulkMolecules
	"""

	geneView = bulkChrmCntr.countsView(kb.geneData['name'])
	geneView.countsInc(1)


def initializeComplexes(bulkMolCntr, kb, randomState, timeStep):
	from wholecell.utils.mc_complexation import mccFormComplexes

	stoichMatrix = kb.complexationStoichMatrix().astype(np.int64, order = "F")

	# Build views

	moleculeNames = kb.complexationMoleculeNames

	molecules = bulkMolCntr.countsView(moleculeNames)

	moleculeCounts = molecules.counts()

	seed = randomState.randint(8**8)

	updatedMoleculeCounts = mccFormComplexes(moleculeCounts, seed, stoichMatrix)

	molecules.countsIs(updatedMoleculeCounts)


def initializeTranscription(bulkMolCntr, uniqueMolCntr, kb, randomState, timeStep):
	"""
	initializeTranscription

	Purpose:
	Initiates transcripts in mid-transcription.

	Method:
	Replaces some RNAP subunits with active RNAP.  These active	RNAPs are in a 
	random position in the elongation process, and are elongating transcripts 
	chosen with respect to their initiation probabilities. Transcript counts 
	are decremented randomly until the RNA mass is approximately its original 
	value.

	Needs attention:
	-Interaction with protein complexes
	-Mass calculations

	"""
	# Calculate the number of possible RNAPs

	inactiveRnap = bulkMolCntr.countView("APORNAP-CPLX[c]")

	activeRnapMax = inactiveRnap.count()

	if activeRnapMax == 0:
		return

	# Calculate the number of RNAPs that should be active

	elngRate = kb.rnaPolymeraseElongationRate.to('nucleotide / s').magnitude

	rnaIds = kb.rnaData["id"]
	rnas = bulkMolCntr.countsView(rnaIds)
	rnaCounts = rnas.counts()
	rnaLengths = kb.rnaData["length"].to("count").magnitude

	rnaLengthAverage = np.dot(rnaCounts, rnaLengths) / rnaCounts.sum()

	fractionActive = 1 - elngRate/rnaLengthAverage * timeStep

	activeRnapCount = np.int64(activeRnapMax * fractionActive)

	# Choose RNAs to polymerize
	rnaCountsPolymerizing = randomState.multinomial(activeRnapCount,
		kb.rnaData["synthProb"])

	# Get the RNA masses

	rnaMasses = (kb.rnaData["mw"].to("fg / mole").magnitude /
		kb.nAvogadro.to("1 / mole").magnitude)

	# Reduce the number of RNAP subunits

	inactiveRnap.countDec(activeRnapCount)

	# Create the lists of RNA indexes

	rnaIndexes = np.empty(activeRnapCount, np.int64)

	startIndex = 0
	for rnaIndex, counts in enumerate(rnaCountsPolymerizing):
		rnaIndexes[startIndex:startIndex+counts] = rnaIndex

		startIndex += counts

	# Choose random RNA lengths

	maxRnaLengths = rnaLengths[rnaIndexes]

	transcriptLengths = (randomState.rand(activeRnapCount) * maxRnaLengths).astype(np.int64)

	# Compute RNA masses
	rnaSequences = kb.transcriptionSequences

	# TODO: standardize this logic w/ process

	ntWeights = np.array([
		345.20, # A
		321.18, # C
		361.20, # G
		322.17, # U
		]) - 17.01 # weight of a hydroxyl

	# TOKB
	hydroxylWeight = 17.01 # counted once for the end of the polymer

	ntWeights *= 1e15/6.022e23
	hydroxylWeight *= 1e15/6.022e23

	transcriptMasses = np.array([
		ntWeights[rnaSequences[rnaIndex, :length]].sum()
		for rnaIndex, length in izip(rnaIndexes, transcriptLengths)
		])

	transcriptMasses[transcriptLengths > 0] += hydroxylWeight

	# Create the unique molecules representations of the ribosomes

	activeRnaPolys = uniqueMolCntr.objectsNew(
		"activeRnaPoly",
		activeRnapCount
		)

	activeRnaPolys.attrIs(
		rnaIndex = rnaIndexes,
		transcriptLength = transcriptLengths,
		massDiff_mRNA = transcriptMasses
		)

	# Compute the amount of RNA mass that needs to be removed
	rnaMassNeeded = transcriptMasses.sum()

	# Using their relative fractions, remove RNAs until the mass is restored

	# TODO: make sure this while-loop isn't too terrible slow

	rnaCountsDecremented = np.zeros_like(rnaCounts)

	while True:
		rnaFractions = rnaCounts / rnaCounts.sum()
		rnaFractionsCumulative = rnaFractions.cumsum()

		rnaIndex = np.where( # sample once from a multinomial distribution
			rnaFractionsCumulative > randomState.rand()
			)[0][0]

		rnaMass = rnaMasses[rnaIndex]

		if rnaMass < rnaMassNeeded:
			rnaMassNeeded -= rnaMass
			rnaCountsDecremented[rnaIndex] += 1
			rnaCounts[rnaIndex] -= 1

		elif rnaMass / 2 < rnaMassNeeded:
			rnaMassNeeded -= rnaMass
			rnaCountsDecremented[rnaIndex] += 1
			rnaCounts[rnaIndex] -= 1

			break

		else:
			break

	rnas.countsDec(rnaCountsDecremented)

def initializeTranslation(bulkMolCntr, uniqueMolCntr, kb, randomState, timeStep):
	"""
	initializeTranslation

	Purpose:
	Initiates peptides in mid-translation.

	Method:
	Replaces some ribosomal subunits with active ribosomes.  These active
	ribosomes are in a random position in the elongation process, and are 
	elongating peptides chosen with respect to the steady-state mRNA
	expression. Protein monomer counts are decremented randomly until the 
	protein mass is approximately its original value.

	Needs attention:
	-Interaction with protein complexes
	-Mass calculations

	"""
	# Calculate the number of possible ribosomes

	subunits = bulkMolCntr.countsView(["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"])

	subunitStoich = np.array([1, 1, 1])

	activeRibosomeMax = (subunits.counts() // subunitStoich).min()

	if activeRibosomeMax == 0:
		return

	# Calculate the number of ribosomes that should be active

	elngRate = kb.ribosomeElongationRate.to('amino_acid / s').magnitude

	monomerIds = kb.monomerData["id"]
	monomers = bulkMolCntr.countsView(monomerIds)
	monomerCounts = monomers.counts()
	monomerLengths = kb.monomerData["length"].to("count").magnitude

	monomerLengthAverage = np.dot(monomerCounts, monomerLengths) / monomerCounts.sum()

	fractionActive = 1 - elngRate/monomerLengthAverage * timeStep

	activeRibosomeCount = np.int64(activeRibosomeMax * fractionActive)

	# Choose proteins to polymerize
	mrnaExpression = kb.rnaExpression["expression"][kb.rnaIndexToMonomerMapping]
	averageMrnaFractions = mrnaExpression/mrnaExpression.sum()

	monomerCountsPolymerizing = randomState.multinomial(activeRibosomeCount, averageMrnaFractions)

	# Compute the current protein mass

	monomerMasses = (kb.monomerData["mw"].to("fg / mole").magnitude /
		kb.nAvogadro.to("1 / mole").magnitude)
	monomerMassTotal = np.dot(monomerCounts, monomerMasses)

	# Reduce the number of ribosome subunits

	subunits.countsDec(activeRibosomeCount * subunitStoich)

	# Create the lists of protein indexes

	proteinIndexes = np.empty(activeRibosomeCount, np.int64)

	startIndex = 0
	for proteinIndex, counts in enumerate(monomerCountsPolymerizing):
		proteinIndexes[startIndex:startIndex+counts] = proteinIndex

		startIndex += counts

	# Choose random protein lengths

	maxPeptideLengths = monomerLengths[proteinIndexes]

	peptideLengths = (randomState.rand(activeRibosomeCount) * maxPeptideLengths).astype(np.int64)

	# Compute peptide masses
	monomerSequences = kb.translationSequences

	aaWeightsIncorporated = kb.translationMonomerWeights

	# TODO: check whether there should be an additional water mass
	peptideMasses = np.array([
		aaWeightsIncorporated[monomerSequences[proteinIndex, :length]].sum()
		for proteinIndex, length in izip(proteinIndexes, peptideLengths)
		])

	# Create the unique molecules representations of the ribosomes

	activeRibosomes = uniqueMolCntr.objectsNew(
		"activeRibosome",
		activeRibosomeCount
		)

	activeRibosomes.attrIs(
		proteinIndex = proteinIndexes,
		peptideLength = peptideLengths,
		massDiff_protein = peptideMasses
		)

	# Compute the amount of monomer mass that needs to be removed
	monomerMassNeeded = peptideMasses.sum()

	# Using their relative fractions, remove monomers until the mass is restored

	# TODO: make sure this while-loop isn't too terrible slow

	monomerCountsDecremented = np.zeros_like(monomerCounts)

	while True:
		monomerFractions = monomerCounts / monomerCounts.sum()
		monomerFractionsCumulative = monomerFractions.cumsum()

		monomerIndex = np.where( # sample once from a multinomial distribution
			monomerFractionsCumulative > randomState.rand()
			)[0][0]

		monomerMass = monomerMasses[monomerIndex]

		if monomerMass < monomerMassNeeded:
			monomerMassNeeded -= monomerMass
			monomerCountsDecremented[monomerIndex] += 1
			monomerCounts[monomerIndex] -= 1

		elif monomerMass / 2 < monomerMassNeeded:
			monomerMassNeeded -= monomerMass
			monomerCountsDecremented[monomerIndex] += 1
			monomerCounts[monomerIndex] -= 1

			break

		else:
			break

	monomers.countsDec(monomerCountsDecremented)

def initializeReplication(uniqueMolCntr, kb):
	'''
	initializeReplication

	Purpose: Create two replication forks represented as unique
	molecules for now at the center of the oriC
	'''
	oricCenter = kb.oriCCenter.to('nucleotide').magnitude
	dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', 4)
	dnaPoly.attrIs(
		chromosomeLocation = np.array([oricCenter, oricCenter, oricCenter, oricCenter]),
		directionIsPositive = np.array([True, True, False, False]),
		isLeading = np.array([True, False, True, False])
		)
