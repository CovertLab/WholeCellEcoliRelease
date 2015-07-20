
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

def initializeBulkChromosome(bulkChrmCntr, kb, randomState, timeStep):
	## Set genes
	initializeGenes(bulkChrmCntr, kb, timeStep)

def initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, kb, randomState, timeStep):
	initializeReplication(uniqueMolCntr, kb)

def initializeProteinMonomers(bulkMolCntr, kb, randomState, timeStep):

	monomersView = bulkMolCntr.countsView(kb.process.translation.monomerData["id"])
	monomerMass = kb.mass.subMass["proteinMass"]
	# TODO: unify this logic with the fitter so it doesn't fall out of step
	# again (look at the calcProteinCounts function)

	monomerExpression = normalize(
		kb.process.transcription.rnaData["expression"][kb.relation.rnaIndexToMonomerMapping] /
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
	rnaMass = kb.mass.subMass["rnaMass"]

	rnaExpression = normalize(kb.process.transcription.rnaData["expression"])

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

	chromosomeView = bulkMolCntr.countsView(kb.moleculeGroups.fullChromosome)
	chromosomeView.countsIs([1])

# TODO: remove checks for zero concentrations (change to assertion)
# TODO: move any rescaling logic to KB/fitting
def initializeSmallMolecules(bulkMolCntr, kb, randomState, timeStep):
	subMass = kb.mass.subMass

	mass = subMass["proteinMass"] + subMass["rnaMass"] + subMass["dnaMass"]

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

def initializeGenes(bulkChrmCntr, kb, timeStep):
	"""
	initializeGenes

	Purpose:
	Initializes the counts of genes in BulkChromosome
	"""

	geneView = bulkChrmCntr.countsView(kb.process.replication.geneData['name'])
	geneView.countsInc(1)

def initializeReplication(uniqueMolCntr, kb):
	'''
	initializeReplication

	Purpose: Create two replication forks represented as unique
	molecules for now at the center of the oriC
	'''
	oricCenter = kb.constants.oriCCenter.asNumber(units.nt)
	dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', 4)
	dnaPoly.attrIs(
		sequenceIdx = np.array([0, 1, 2, 3]),
		sequenceLength = np.array([0, 0, 0, 0]),
		replicationRound = np.array([0, 0, 0, 0]),
		replicationDivision = np.array([0, 0, 0, 0])
		)
	# TODO: DELETE!
	# dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', 8)
	# dnaPoly.attrIs(
	# 	sequenceIdx = np.array([0, 1, 2, 3, 0, 1, 2, 3]),
	# 	#sequenceLength = np.array([0, 0, 0, 0]),
	# 	sequenceLength = np.array([229500.0, 229500.0, 229500.0, 229500.0, 229500.0, 229500.0, 229500.0, 229500.0]),
	# 	replicationRound = np.array([1, 1, 1, 1, 1, 1, 1, 1]),
	# 	replicationDivision = np.array([0, 0, 0, 0, 1, 1, 1, 1])
	# 	)

	# dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', 16)
	# dnaPoly.attrIs(
	# 	sequenceIdx = np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]),
	# 	#sequenceLength = np.array([0, 0, 0, 0]),
	# 	sequenceLength = np.array([22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0, 22950.0]),
	# 	replicationRound = np.array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]),
	# 	replicationDivision = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1])
	# 	)

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
