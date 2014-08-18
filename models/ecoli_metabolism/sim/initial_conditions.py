
from __future__ import division

from itertools import izip

import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.fitter import countsFromMassAndExpression
from reconstruction.ecoli.fitter import normalize
from reconstruction.ecoli.fitter import calcProteinCounts
from wholecell.utils import units

def calcInitialConditions(sim, kb):
	randomState = sim.randomState

	timeStep = sim.timeStepSec() # This is a poor solution but will suffice for now

	bulkMolCntr = sim.states['BulkMolecules'].container

	# Set up states
	initializeBulkMolecules(bulkMolCntr, kb, randomState, timeStep)


def initializeBulkMolecules(bulkMolCntr, kb, randomState, timeStep):

	# Initialize protein
	initializeProtein(bulkMolCntr, kb, randomState, timeStep)

	# Initialize RNA
	initializeRNA(bulkMolCntr, kb, randomState, timeStep)

	# Initialize DNA
	initializeDNA(bulkMolCntr, kb, randomState, timeStep)

	# Set other biomass components
	initializeBulkComponents(bulkMolCntr, kb, randomState, timeStep)


def initializeProtein(bulkMolCntr, kb, randomState, timeStep):
	polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedAA_IDs]

	proteinComposition = kb.monomerData["aaCounts"].asNumber()

	initialDryMass = kb.avgCellDryMassInit

	proteinMassFraction = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].asUnit(units.min).asNumber() == 60.0
		]["proteinMassFraction"]

	initialProteinMass = initialDryMass * proteinMassFraction

	proteinCounts = calcProteinCounts(kb, initialProteinMass)

	polymerizedCounts = np.dot(proteinComposition.T, proteinCounts)
	
	bulkMolCntr.countsIs(
		np.int64(polymerizedCounts),
		polymerizedIDs
		)


def initializeRNA(bulkMolCntr, kb, randomState, timeStep):
	# TODO: move duplicate logic to the KB

	polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedNT_IDs]

	## Find the average RNA composition

	synthProb = kb.rnaData["synthProb"]
	compositionAll = kb.rnaData["countsACGU"].asNumber()

	# TODO: better model the variance of this distribution

	compositionUnnormed = np.dot(compositionAll.T, synthProb)

	monomerComposition = compositionUnnormed / compositionUnnormed.sum()

	## Find the average total transcription rate with respect to cell age

	initialDryMass = kb.avgCellDryMassInit.asUnit(units.fg).asNumber()

	rnaMassFraction = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].asUnit(units.min).asNumber() == 60.0
		]["rnaMassFraction"]

	initialRnaMass = initialDryMass * rnaMassFraction

	monomerMWs = kb.transcriptionMonomerWeights

	initialMonomerCounts = np.int64(
		initialRnaMass * monomerComposition / monomerMWs
		)
	
	bulkMolCntr.countsIs(
		initialMonomerCounts,
		polymerizedIDs
		)


def initializeDNA(bulkMolCntr, kb, randomState, timeStep):
	# TODO: move duplicate logic to the KB
	
	polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedDNT_IDs]

	bulkMolCntr.countsIs(
		[
			kb.genomeSeq.count("A") + kb.genomeSeq.count("T"),
			kb.genomeSeq.count("C") + kb.genomeSeq.count("G"),
			kb.genomeSeq.count("G") + kb.genomeSeq.count("C"),
			kb.genomeSeq.count("T") + kb.genomeSeq.count("A")
		],
		polymerizedIDs
		)


def initializeBulkComponents(bulkMolCntr, kb, randomState, timeStep):

	massFractions = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].asUnit(units.min).asNumber() == 60.0
		].fullArray()

	initDryMass = kb.avgCellDryMassInit.asUnit(units.g).asNumber()

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

	bulkMolCntr.countsIs(
		countsToAdd,
		poolIds
		)
