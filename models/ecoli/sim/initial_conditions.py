
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
import os

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.utils.fitting import normalize, countsFromMassAndExpression, calcProteinCounts
from wholecell.utils import units

from wholecell.io.tablereader import TableReader

def calcInitialConditions(sim, kb):
	assert sim._inheritedStatePath == None
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
	initializeSmallMolecules(bulkMolCntr, kb, randomState, timeStep)

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

	monomersView = bulkMolCntr.countsView(kb.process.translation.monomerData["id"])
	monomerMass = kb.mass.massFractions["proteinMass"]
	# TODO: unify this logic with the fitter so it doesn't fall out of step
	# again (look at the calcProteinCounts function)

	monomerExpression = normalize(
		kb.process.transcription.rnaData['expression'][kb.relation.rnaIndexToMonomerMapping] /
		(np.log(2) / kb.doubling_time.asNumber(units.s) + kb.process.translation.monomerData["degRate"].asNumber(1 / units.s))
		)

	nMonomers = countsFromMassAndExpression(
		monomerMass.asNumber(units.g),
		kb.process.translation.monomerData["mw"].asNumber(units.g/units.mol),
		monomerExpression,
		kb.constants.nAvogadro.asNumber(1/units.mol)
		)

	monomersView.countsIs(
		randomState.multinomial(nMonomers, monomerExpression)
		)

def initializeRNA(bulkMolCntr, kb, randomState, timeStep):

	rnaView = bulkMolCntr.countsView(kb.process.transcription.rnaData["id"])
	rnaMass = kb.mass.massFractions["rnaMass"]

	rnaExpression = normalize(kb.process.transcription.rnaData['expression'])

	nRnas = countsFromMassAndExpression(
		rnaMass.asNumber(units.g),
		kb.process.transcription.rnaData["mw"].asNumber(units.g/units.mol),
		rnaExpression,
		kb.constants.nAvogadro.asNumber(1/units.mol)
		)

	rnaView.countsIs(
		randomState.multinomial(nRnas, rnaExpression)
		)

def initializeDNA(bulkMolCntr, kb, randomState, timeStep):

	polymerizedView = bulkMolCntr.countsView([id_ + "[c]" for id_ in kb.moleculeGroups.polymerizedDNT_IDs])

	polymerizedView.countsIs([
		kb.process.replication.genome_A_count + kb.process.replication.genome_T_count,
		kb.process.replication.genome_C_count + kb.process.replication.genome_G_count,
		kb.process.replication.genome_G_count + kb.process.replication.genome_C_count,
		kb.process.replication.genome_T_count + kb.process.replication.genome_A_count
		])

def initializeSmallMolecules(bulkMolCntr, kb, randomState, timeStep):
	massFractions60 = kb.mass.massFractions

	mass = massFractions60["proteinMass"] + massFractions60["rnaMass"] + massFractions60["dnaMass"]

	# We have to remove things with zero concentration because taking the inverse of zero isn't so nice.
	poolIds = [x for idx, x in enumerate(kb.process.metabolism.metabolitePoolIDs) if kb.process.metabolism.metabolitePoolConcentrations.asNumber()[idx] > 0]
	poolConcentrations = (units.mol / units.L) * np.array([x for x in kb.process.metabolism.metabolitePoolConcentrations.asNumber() if x > 0])

	cellDensity = kb.constants.cellDensity
	mws = kb.getter.getMass(poolIds)
	concentrations = poolConcentrations.copy()

	diag = (cellDensity / (mws * concentrations) - 1).asNumber()
	A = -1 * np.ones((diag.size, diag.size))
	A[np.diag_indices(diag.size)] = diag
	b = mass.asNumber(units.g) * np.ones(diag.size)

	massesToAdd = units.g * np.linalg.solve(A, b)
	countsToAdd = massesToAdd / mws * kb.constants.nAvogadro

	V = (mass + units.sum(massesToAdd)) / cellDensity

	assert np.allclose(
		(countsToAdd / kb.constants.nAvogadro / V).asNumber(units.mol / units.L),
		(poolConcentrations).asNumber(units.mol / units.L)
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		poolIds
		)

	# GDP POOL
	ribosomeSubunits = bulkMolCntr.countsView(
		np.hstack(
			(kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitIds'], kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitIds'])
			)
		)
	ribosomeSubunitStoich = np.hstack(
			(kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitStoich'], kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitStoich'])
			)

	activeRibosomeMax = (ribosomeSubunits.counts() // ribosomeSubunitStoich).min()
	elngRate = kb.constants.ribosomeElongationRate.asNumber(units.aa / units.s)
	T_d = kb.doubling_time.asNumber(units.s)
	dt = kb.constants.timeStep.asNumber(units.s)

	activeRibosomesLastTimeStep = activeRibosomeMax * np.exp( np.log(2) / T_d * (T_d - dt)) / 2
	gtpsHydrolyzedLastTimeStep = activeRibosomesLastTimeStep * elngRate * kb.constants.gtpPerTranslation

	bulkMolCntr.countsInc(
		gtpsHydrolyzedLastTimeStep,
		["GDP[c]"]
		)

def initializeGenes(bulkChrmCntr, kb, timeStep):
	"""
	initializeGenes

	Purpose:
	Initializes the counts of genes in BulkChromosome
	"""

	geneView = bulkChrmCntr.countsView(kb.process.replication.geneData['name'])
	geneView.countsInc(1)

def initializeComplexes(bulkMolCntr, kb, randomState, timeStep):
	from wholecell.utils.mc_complexation import mccFormComplexes

	stoichMatrix = kb.process.complexation.stoichMatrix().astype(np.int64, order = "F")

	# Build views

	moleculeNames = kb.process.complexation.moleculeNames

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

	activeRnapCount = np.int64(activeRnapMax * kb.fracActiveRnap * 2) # HACK

	# Load RNA data
	rnaIds = kb.process.transcription.rnaData["id"]
	rnas = bulkMolCntr.countsView(rnaIds)
	rnaCounts = rnas.counts()
	rnaLengths = kb.process.transcription.rnaData["length"].asNumber(units.count)

	rnaMasses = (kb.process.transcription.rnaData["mw"].asNumber(units.fg / units.mol) /
		kb.constants.nAvogadro.asNumber(1 / units.mol))

	# Choose RNAs to polymerize
	rnaCountsPolymerizing = randomState.multinomial(activeRnapCount,
		kb.process.transcription.rnaData["synthProb"])

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
	rnaSequences = kb.process.transcription.transcriptionSequences

	# TODO: standardize this logic w/ process

	ntWeights = kb.process.transcription.transcriptionMonomerWeights

	# TOKB

	transcriptMasses = np.array([
		ntWeights[rnaSequences[rnaIndex, :length]].sum()
		for rnaIndex, length in izip(rnaIndexes, transcriptLengths)
		])

	transcriptMasses[transcriptLengths > 0] += kb.process.transcription.transcriptionEndWeight

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

	ribosomeSubunits = bulkMolCntr.countsView([kb.moleculeGroups.s30_fullComplex[0], kb.moleculeGroups.s50_fullComplex[0]])
	ribosomeSubunitStoich = np.array([1,1])
	activeRibosomeMax = (ribosomeSubunits.counts() // ribosomeSubunitStoich).min()

	if activeRibosomeMax == 0:
		return

	# Calculate the number of ribosomes that should be active

	elngRate = kb.constants.ribosomeElongationRate.asNumber(units.aa / units.s)

	monomerIds = kb.process.translation.monomerData["id"]
	monomers = bulkMolCntr.countsView(monomerIds)
	monomerCounts = monomers.counts()
	monomerLengths = kb.process.translation.monomerData["length"].asNumber(units.count)

	monomerLengthAverage = np.dot(monomerCounts, monomerLengths) / monomerCounts.sum()

	fractionActive = 1 - elngRate/monomerLengthAverage * timeStep

	activeRibosomeCount = np.int64(activeRibosomeMax * fractionActive)

	# Choose proteins to polymerize
	mrnaExpression = kb.process.transcription.rnaData["expression"][kb.relation.rnaIndexToMonomerMapping]
	averageMrnaFractions = mrnaExpression/mrnaExpression.sum()

	monomerCountsPolymerizing = randomState.multinomial(activeRibosomeCount, averageMrnaFractions)

	# Compute the current protein mass

	monomerMasses = (kb.process.translation.monomerData["mw"].asNumber(units.fg / units.mol) /
		kb.constants.nAvogadro.asNumber(1 / units.mol))
	monomerMassTotal = np.dot(monomerCounts, monomerMasses)

	# Reduce the number of ribosome subunits

	ribosomeSubunits.countsDec(activeRibosomeCount * ribosomeSubunitStoich)

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
	monomerSequences = kb.process.translation.translationSequences

	aaWeightsIncorporated = kb.process.translation.translationMonomerWeights

	peptideMasses = np.array([
		aaWeightsIncorporated[monomerSequences[proteinIndex, :length]].sum()
		for proteinIndex, length in izip(proteinIndexes, peptideLengths)
		])

	peptideMasses[peptideLengths > 0] += kb.process.translation.translationEndWeight

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
	oricCenter = kb.constants.oriCCenter.asNumber(units.nt)
	dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', 4)
	dnaPoly.attrIs(
		chromosomeLocation = np.array([oricCenter, oricCenter, oricCenter, oricCenter]),
		directionIsPositive = np.array([True, True, False, False]),
		isLeading = np.array([True, False, True, False])
		)

def setDaughterInitialConditions(sim, kb):
	assert sim._inheritedStatePath != None

	bulk_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "BulkMolecules"))
	sim.states["BulkMolecules"].tableLoad(bulk_table_reader, 0)

	bulk_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "BulkChromosome"))
	sim.states["BulkChromosome"].tableLoad(bulk_table_reader, 0)

	unique_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "UniqueMolecules"))
	sim.states["UniqueMolecules"].tableLoad(unique_table_reader, 0)

	time_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "Time"))
	initialTime = TableReader(os.path.join(sim._inheritedStatePath, "Time")).readAttribute("initialTime")
	sim._initialTime = initialTime

	# TODO: This is a hack until we actually get a replication initialization process working
	if len(sim.states['UniqueMolecules'].container.objectsInCollection('dnaPolymerase')) == 0:
		initializeReplication(sim.states["UniqueMolecules"].container, kb)

