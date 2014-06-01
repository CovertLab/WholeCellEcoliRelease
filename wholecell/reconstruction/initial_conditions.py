
from __future__ import division

from itertools import izip

import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.reconstruction.fitter import countsFromMassAndExpression
from wholecell.reconstruction.fitter import normalize

def calcInitialConditions(sim, kb):
	randStream = sim.randStream

	bulk = sim.states['BulkMolecules']
	unique = sim.states["UniqueMolecules"]

	# Set up states
	initializeBulk(bulk.container, kb, randStream)

	# Modify states for specific processes
	# initializeTranscription(bulk.container, unique.container, kb, randStream)
	initializeTranslation(bulk.container, unique.container, kb, randStream)


def initializeBulk(bulkContainer, kb, randStream):

	## Set protein counts from expression
	initializeProteinMonomers(bulkContainer, kb, randStream)

	## Set RNA counts from expression
	initializeRNA(bulkContainer, kb, randStream)

	## Set DNA
	initializeDNA(bulkContainer, kb, randStream)

	## Set other biomass components
	initializeBulkComponents(bulkContainer, kb, randStream)

	## Set pools
	initializePools(bulkContainer, kb, randStream)

	## Set water
	initializeBulkWater(kb, bulkContainer, randStream)


def initializeProteinMonomers(bulkContainer, kb, randStream):
	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].to('min').magnitude == 60
		]

	monomersView = bulkContainer.countsView(kb.monomerData["id"])
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
		randStream.mnrnd(nMonomers, monomerExpression)
		)

	# monomersView.countsIs(nMonomers * monomerExpression)


def initializeRNA(bulkContainer, kb, randStream):
	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].magnitude == 60
		]

	rnaView = bulkContainer.countsView(kb.rnaData["id"])
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
		randStream.mnrnd(nRnas, rnaExpression)
		)

	# rnaView.countsIs(nRnas * rnaExpression)


def initializeDNA(bulkContainer, kb, randStream):

	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].magnitude == 60
		]

	dNTPs = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
	dnmpIds = ["DAMP[n]", "DCMP[n]", "DGMP[n]", "DTMP[n]"]

	dnmpsView = bulkContainer.countsView(dnmpIds)
	dnaMassFraction = float(dryComposition60min["dnaMassFraction"])
	dnaMass = kb.avgCellDryMassInit.magnitude * dnaMassFraction

	dnaExpression = normalize(np.array([
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		], dtype = np.float64))

	mws = np.array([
		kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == x][0].magnitude for x in dnmpIds]
		) # This is a hack. Without a real chromosome, though, it's all a hack

	nDntps = countsFromMassAndExpression(
		dnaMass,
		mws,
		dnaExpression,
		kb.nAvogadro.magnitude
		)

	dnmpsView.countsIs(
		randStream.mnrnd(nDntps, dnaExpression)
		)


def initializeBulkComponents(bulkContainer, kb, randStream):
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

	notPRDBulkView = bulkContainer.countsView(notPRDMetabolites)

	notPRDBiomassView = biomassContainer.countsView(notPRDMetabolites)

	notPRDBulkView.countsIs((
		kb.avgCellDryMassInit.to("DCW_gram").magnitude *
		notPRDBiomassView.counts() *
		kb.nAvogadro.to("1 / millimole").magnitude
		))


def initializePools(bulkContainer, kb, randStream):
	# Note: This is adding dry biomass, so the cell will appear heavier

	from wholecell.reconstruction.knowledge_base_ecoli import AMINO_ACID_1_TO_3_ORDERED

	biomassContainer = BulkObjectsContainer(
		list(kb.wildtypeBiomass["metaboliteId"]), dtype = np.dtype("float64")
		)
	biomassContainer.countsIs(
		kb.wildtypeBiomass["biomassFlux"].to("millimole/DCW_gram").magnitude
		)

	ntpIds = ["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]
	dntpIds = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]
	aaIds = [x for x in AMINO_ACID_1_TO_3_ORDERED.values() if x != "SEC-L[c]"]

	ntpsBiomassView = biomassContainer.countsView(ntpIds)
	dntpsBiomassView = biomassContainer.countsView(dntpIds)
	aasBiomassView = biomassContainer.countsView(aaIds)

	ppiBulkView = bulkContainer.countView("PPI[c]")
	ntpsBulkView = bulkContainer.countsView(ntpIds)
	dntpsBulkView = bulkContainer.countsView(dntpIds)
	aasBulkView = bulkContainer.countsView(aaIds)

	dt = kb.timeStep.to("second").magnitude
	tau_d = kb.cellCycleLen.to("second").magnitude

	## NTPs
	ntpsFromLastStep = (
		ntpsBiomassView.counts() *
		(1 - np.exp(-np.log(2) / tau_d * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)
	ntpsBulkView.countsIs(ntpsFromLastStep)

	## dNTPs
	dntpsFromLastStep = (
		dntpsBiomassView.counts() * (1 - np.exp(-np.log(2) / tau_d * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)
	dntpsBulkView.countsIs(dntpsFromLastStep)

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

	dntpsFromTwoStepsAgo = (
		dntpsBiomassView.counts() *
		(np.exp(-np.log(2) / tau_d * dt) - np.exp(-np.log(2) / tau_d * 2 * dt)) *
		kb.nAvogadro.to("1 / millimole").magnitude *
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		)

	ppiFromNtps = np.round(np.sum(
		ntpsFromLastStep
		))

	ppiFromDntps = np.round(np.sum(
		dntpsFromLastStep
		))

	ppiBulkView.countIs(ppiFromNtps + ppiFromDntps)


def initializeBulkWater(kb, bulkContainer, randStream):
	h2oView = bulkContainer.countView('H2O[c]')

	nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
	mwH2O = kb.bulkMolecules["mass"][kb.bulkMolecules["moleculeId"] == "H2O[c]"].magnitude
	avgCellWaterMassInit = kb.avgCellWaterMassInit.to('water_g').magnitude

	h2oView.countIs(
		(avgCellWaterMassInit) / mwH2O * nAvogadro
		)

def initializeTranslation(bulkContainer, uniqueContainer, kb, randStream):
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

	subunits = bulkContainer.countsView(["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"])

	subunitStoich = np.array([1, 1, 1])

	activeRibosomeMax = (subunits.counts() * subunitStoich).min()

	# Calculate the number of ribosomes that should be active

	elngRate = kb.ribosomeElongationRate.to('amino_acid / s').magnitude

	monomerIds = kb.monomerData["id"]
	monomers = bulkContainer.countsView(monomerIds)
	monomerCounts = monomers.counts()
	monomerLengths = kb.monomerData["length"].to("count").magnitude

	monomerLengthAverage = np.dot(monomerCounts, monomerLengths) / monomerCounts.sum()

	fractionActive = 1 - elngRate/monomerLengthAverage

	activeRibosomeCount = np.int64(activeRibosomeMax * fractionActive)

	# Choose proteins to polymerize
	mrnaExpression = kb.rnaExpression["expression"][kb.rnaIndexToMonomerMapping]
	averageMrnaFractions = mrnaExpression/mrnaExpression.sum()

	monomerCountsPolymerizing = randStream.mnrnd(activeRibosomeCount, averageMrnaFractions)

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

	peptideLengths = (randStream.rand(activeRibosomeCount) * maxPeptideLengths).astype(np.int64)

	# Compute peptide masses
	monomerSequences = kb.translationSequences

	# TODO: standardize this logic w/ process

	h2oWeight = (
		kb.bulkMolecules[
			kb.bulkMolecules["moleculeId"] == "H2O[c]"
			]["mass"].to("fg / mole").magnitude /
		kb.nAvogadro.to("1 / mole").magnitude
		)

	aaWeights = np.array([
		kb.bulkMolecules[
			kb.bulkMolecules["moleculeId"] == x
			]["mass"].to("fg / mole").magnitude /
		kb.nAvogadro.to("1 / mole").magnitude
		for x in kb.aaIDs
		if len(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == x]["mass"])
		])

	aaWeightsIncorporated = aaWeights - h2oWeight

	# TODO: check whether there should be an additional water mass
	peptideMasses = np.array([
		aaWeightsIncorporated[monomerSequences[proteinIndex, :length]].sum()
		for proteinIndex, length in izip(proteinIndexes, peptideLengths)
		])

	# Create the unique molecules representations of the ribosomes

	activeRibosomes = uniqueContainer.objectsNew(
		"activeRibosome",
		activeRibosomeCount
		)

	activeRibosomes.attrIs(
		proteinIndex = proteinIndexes,
		peptideLength = peptideLengths,
		massDiffProtein = peptideMasses
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
			monomerFractionsCumulative > randStream.rand()
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
