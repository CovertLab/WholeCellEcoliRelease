
from __future__ import division

from itertools import izip

import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.utils.fitting import normalize, countsFromMassAndExpression
from reconstruction.ecoli.fitter import calcProteinCounts
from reconstruction.ecoli.compendium import growth_data
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
	initializeSmallMolecules(bulkMolCntr, kb, randomState, timeStep)


def initializeProtein(bulkMolCntr, kb, randomState, timeStep):
	polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedAA_IDs]

	proteinComposition = kb.monomerData["aaCounts"].asNumber()

	g = growth_data.GrowthData(kb)
	initialProteinMass = g.massFractions(60)["proteinMass"]

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

	g = growth_data.GrowthData(kb)
	initialRnaMass = g.massFractions(60)["rnaMass"].asNumber(units.fg)

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
			kb.genome_A_count + kb.genome_T_count,
			kb.genome_C_count + kb.genome_G_count,
			kb.genome_G_count + kb.genome_C_count,
			kb.genome_T_count + kb.genome_A_count
		],
		polymerizedIDs
		)


def initializeSmallMolecules(bulkMolCntr, kb, randomState, timeStep):

	g = growth_data.GrowthData(kb)
	massFractions60 = g.massFractions(60)

	mass = massFractions60["proteinMass"] + massFractions60["rnaMass"] + massFractions60["dnaMass"]

	# We have to remove things with zero concentration because taking the inverse of zero isn't so nice.
	poolIds = [x for idx, x in enumerate(kb.metabolitePoolIDs) if kb.metabolitePoolConcentrations.asNumber()[idx] > 0]
	poolConcentrations = (units.mol / units.L) * np.array([x for x in kb.metabolitePoolConcentrations.asNumber() if x > 0])

	cellDensity = kb.cellDensity
	mws = kb.getMass(poolIds)
	concentrations = poolConcentrations.copy()

	diag = (cellDensity / (mws * concentrations) - 1).asNumber()
	A = -1 * np.ones((diag.size, diag.size))
	A[np.diag_indices(diag.size)] = diag
	b = mass.asNumber(units.g) * np.ones(diag.size)

	massesToAdd = units.g * np.linalg.solve(A, b)
	countsToAdd = massesToAdd / mws * kb.nAvogadro

	V = (mass + units.sum(massesToAdd)) / cellDensity

	assert np.allclose(
		(countsToAdd / kb.nAvogadro / V).asNumber(units.mol / units.L),
		(poolConcentrations).asNumber(units.mol / units.L)
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		poolIds
		)