
"""

TODO:
- document math
- replace fake metabolite pools with measured metabolite pools
- raise/warn if physiological metabolite pools appear to be smaller than what
 is needed at this time step size

"""

from __future__ import division

from itertools import izip

import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.fitter import countsFromMassAndExpression
from reconstruction.ecoli.fitter import normalize
from wholecell.utils import units

def calcInitialConditions(sim, kb):
	randomState = sim.randomState

	timeStep = sim.timeStepSec() # This is a poor solution but will suffice for now

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
		kb.cellDryMassComposition["doublingTime"].asUnit(units.min).asNumber() == 60
		]

	monomersView = bulkMolCntr.countsView(kb.monomerData["id"])
	monomerMassFraction = float(dryComposition60min["proteinMassFraction"])
	monomerMass = kb.avgCellDryMassInit.asUnit(units.g) * monomerMassFraction

	# TODO: unify this logic with the fitter so it doesn't fall out of step 
	# again (look at the calcProteinCounts function)

	monomerExpression = normalize(
		kb.rnaExpression['expression'][kb.rnaIndexToMonomerMapping] /
		(np.log(2) / kb.cellCycleLen.asUnit(units.s).asNumber() + kb.monomerData["degRate"].asUnit(1 / units.s).asNumber())
		)

	nMonomers = countsFromMassAndExpression(
		monomerMass.asUnit(units.g).asNumber(),
		kb.monomerData["mw"].asUnit(units.g/units.mol).asNumber(),
		monomerExpression,
		kb.nAvogadro.asUnit(1/units.mol).asNumber()
		)

	monomersView.countsIs(
		randomState.multinomial(nMonomers, monomerExpression)
		)

	# monomersView.countsIs(nMonomers * monomerExpression)

	## Uncomment for debugging
	# M = kb.getMass([monomersView._container._objectNames[x] for x in monomersView._indexes]).asNumber()
	# initDryMass = kb.avgCellDryMassInit.to("DCW_gram").asNumber()
	# print "Expected Protein Mass: %g" % (initDryMass * dryComposition60min.fullArray()["proteinMassFraction"])[0]
	# print "Actual Protein Mass: %g" % (monomersView.counts() / kb.nAvogadro.asNumber() * M).sum()


def initializeRNA(bulkMolCntr, kb, randomState, timeStep):
	dryComposition60min = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].asNumber() == 60
		]

	rnaView = bulkMolCntr.countsView(kb.rnaData["id"])
	rnaMassFraction = float(dryComposition60min["rnaMassFraction"])
	rnaMass = kb.avgCellDryMassInit * rnaMassFraction

	rnaExpression = normalize(kb.rnaExpression['expression'])

	nRnas = countsFromMassAndExpression(
		rnaMass.asUnit(units.g).asNumber(),
		kb.rnaData["mw"].asUnit(units.g/units.mol).asNumber(),
		rnaExpression,
		kb.nAvogadro.asUnit(1/units.mol).asNumber()
		)

	rnaView.countsIs(
		randomState.multinomial(nRnas, rnaExpression)
		)

	# rnaView.countsIs(nRnas * rnaExpression)

	## Uncomment for debugging
	# M = kb.getMass([rnaView._container._objectNames[x] for x in rnaView._indexes]).asNumber()
	# initDryMass = kb.avgCellDryMassInit.to("DCW_gram").asNumber()
	# print "Expected RNA Mass: %g" % (initDryMass * dryComposition60min.fullArray()["rnaMassFraction"])[0]
	# print "Actual RNA Mass: %g" % (rnaView.counts() / kb.nAvogadro.asNumber() * M).sum()



def initializeDNA(bulkMolCntr, kb, randomState, timeStep):

	polymerizedView = bulkMolCntr.countsView([id_ + "[c]" for id_ in kb.polymerizedDNT_IDs])

	polymerizedView.countsIs([
		kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
		kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
		kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
		kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		])


	## Uncomment for debugging
	# M = kb.getMass([dnmpsView._container._objectNames[x] for x in dnmpsView._indexes]).asNumber()
	# initDryMass = kb.avgCellDryMassInit.to("DCW_gram").asNumber()
	# print "Expected DNA Mass: %g" % (initDryMass * dryComposition60min.fullArray()["dnaMassFraction"])[0]
	# print "Actual DNA Mass: %g" % (dnmpsView.counts() / kb.nAvogadro.asNumber() * M).sum()

def initializeBulkComponents(bulkMolCntr, kb, randomState, timeStep):

	massFractions = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].asUnit(units.min).asNumber() == 60.0
		].fullArray()

	initDryMass = kb.avgCellDryMassInit.asUnit(units.g).asNumber()
	cellMass = (
		kb.avgCellDryMassInit.asUnit(units.g).asNumber()
		# + kb.avgCellWaterMassInit.asNumber()
		)

	poolIds = kb.metabolitePoolIDs[:]

	mass = initDryMass
	mass -= massFractions["glycogenMassFraction"] * initDryMass
	mass -= massFractions["mureinMassFraction"] * initDryMass
	mass -= massFractions["lpsMassFraction"] * initDryMass
	mass -= massFractions["lipidMassFraction"] * initDryMass
	mass -= massFractions["inorganicIonMassFraction"] * initDryMass
	mass -= massFractions["solublePoolMassFraction"] * initDryMass

	# We have to remove things with zero concentration because taking the inverse of zero isn't so nice.
	poolIds = [x for idx, x in enumerate(kb.metabolitePoolIDs) if kb.metabolitePoolConcentrations.asNumber()[idx] > 0]
	poolConcentrations = np.array([x for x in kb.metabolitePoolConcentrations.asNumber() if x > 0])

	cellVolume = cellMass / kb.cellDensity
	cellDensity = kb.cellDensity.asUnit(units.g / units.L).asNumber()
	mws = kb.getMass(poolIds).asUnit(units.g / units.mol).asNumber()
	concentrations = poolConcentrations.copy()

	diag = cellDensity / (mws * concentrations) - 1
	A = -1 * np.ones((diag.size, diag.size))
	A[np.diag_indices(diag.size)] = diag
	b = mass * np.ones(diag.size)


	massesToAdd = np.linalg.solve(A, b)
	countsToAdd = massesToAdd / mws * kb.nAvogadro.asUnit(1 / units.mol).asNumber()

	V = (mass + massesToAdd.sum()) / cellDensity

	assert np.allclose(countsToAdd / kb.nAvogadro.asNumber() / V, poolConcentrations)

	## Uncomment the following if you want to look at how well our two different accountings of mass agree
	# M = kb.getMass(bulkMolCntr._objectNames)
	# (bulkMolCntr.counts() / kb.nAvogadro.asNumber() * M).sum()
	# (bulkMolCntr.counts() / kb.nAvogadro.asNumber() * M).sum() / (mass)
	# import ipdb; ipdb.set_trace()

	bulkMolCntr.countsIs(
		countsToAdd,
		poolIds
		)

	## Uncomment the following if you want to look at how well our two different accountings of mass agree
	# M = kb.getMass(bulkMolCntr._objectNames)
	# (bulkMolCntr.counts() / kb.nAvogadro.asNumber() * M).sum()
	# (bulkMolCntr.counts() / kb.nAvogadro.asNumber() * M).sum() / (mass + massesToAdd.sum())
	# import ipdb; ipdb.set_trace()



	subunits = bulkMolCntr.countsView(["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"])
	subunitStoich = np.array([1, 1, 1])
	activeRibosomeMax = (subunits.counts() // subunitStoich).min()
	elngRate = kb.ribosomeElongationRate.asUnit(units.aa / units.s).asNumber()
	T_d = kb.cellCycleLen.asUnit(units.s).asNumber()
	dt = kb.timeStep.asUnit(units.s).asNumber()

	activeRibosomesLastTimeStep = activeRibosomeMax * np.exp( np.log(2) / T_d * (T_d - dt)) / 2
	gtpsHydrolyzedLastTimeStep = activeRibosomesLastTimeStep * elngRate * kb.gtpPerTranslation

	bulkMolCntr.countsInc(
		gtpsHydrolyzedLastTimeStep,
		["GDP[c]"]
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

	elngRate = kb.rnaPolymeraseElongationRate.asUnit(units.nt / units.s).asNumber()

	rnaIds = kb.rnaData["id"]
	rnas = bulkMolCntr.countsView(rnaIds)
	rnaCounts = rnas.counts()
	rnaLengths = kb.rnaData["length"].asUnit(units.count).asNumber()

	rnaLengthAverage = np.dot(rnaCounts, rnaLengths) / rnaCounts.sum()

	fractionActive = 1 - elngRate/rnaLengthAverage * timeStep

	activeRnapCount = np.int64(activeRnapMax * fractionActive)

	# Choose RNAs to polymerize
	rnaCountsPolymerizing = randomState.multinomial(activeRnapCount,
		kb.rnaData["synthProb"])

	# Get the RNA masses

	rnaMasses = (kb.rnaData["mw"].asUnit(units.fg / units.mol).asNumber() /
		kb.nAvogadro.asUnit(1 / units.mol).asNumber())

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

	ntWeights = kb.transcriptionMonomerWeights

	# TOKB

	transcriptMasses = np.array([
		ntWeights[rnaSequences[rnaIndex, :length]].sum()
		for rnaIndex, length in izip(rnaIndexes, transcriptLengths)
		])

	transcriptMasses[transcriptLengths > 0] += kb.transcriptionEndWeight

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

	elngRate = kb.ribosomeElongationRate.asUnit(units.aa / units.s).asNumber()

	monomerIds = kb.monomerData["id"]
	monomers = bulkMolCntr.countsView(monomerIds)
	monomerCounts = monomers.counts()
	monomerLengths = kb.monomerData["length"].asUnit(units.count).asNumber()

	monomerLengthAverage = np.dot(monomerCounts, monomerLengths) / monomerCounts.sum()

	fractionActive = 1 - elngRate/monomerLengthAverage * timeStep

	activeRibosomeCount = np.int64(activeRibosomeMax * fractionActive)

	# Choose proteins to polymerize
	mrnaExpression = kb.rnaExpression["expression"][kb.rnaIndexToMonomerMapping]
	averageMrnaFractions = mrnaExpression/mrnaExpression.sum()

	monomerCountsPolymerizing = randomState.multinomial(activeRibosomeCount, averageMrnaFractions)

	# Compute the current protein mass

	monomerMasses = (kb.monomerData["mw"].asUnit(units.fg / units.mol).asNumber() /
		kb.nAvogadro.asUnit(1 / units.mol).asNumber())
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

	peptideMasses = np.array([
		aaWeightsIncorporated[monomerSequences[proteinIndex, :length]].sum()
		for proteinIndex, length in izip(proteinIndexes, peptideLengths)
		])

	peptideMasses[peptideLengths > 0] += kb.translationEndWeight

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
	oricCenter = kb.oriCCenter.asUnit(units.nt).asNumber()
	dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', 4)
	dnaPoly.attrIs(
		chromosomeLocation = np.array([oricCenter, oricCenter, oricCenter, oricCenter]),
		directionIsPositive = np.array([True, True, False, False]),
		isLeading = np.array([True, False, True, False])
		)
