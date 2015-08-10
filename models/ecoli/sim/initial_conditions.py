
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
from wholecell.utils.polymerize import buildSequences, computeMassIncrease
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
	"""
	initializeReplication

	Purpose: Create the appropriate number of replication forks given the cell growth rate.
	"""

	## Determine the number and location of replication forks at the start of the cell cycle
	# Find growth rate constants
	C = kb.constants.c_period.asUnit(units.min)
	D = kb.constants.d_period.asUnit(units.min)
	tau = kb.doubling_time.asUnit(units.min)
	genome_length = kb.process.replication.genome_length
	replication_length = np.ceil(.5*genome_length) * units.nt

	# Generate arrays specifying appropriate replication conditions
	sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

	# Determine the number of OriC's currently present in the cell
	numOric = determineNumOriC(C, D, tau)

	# Return if no replication is occuring at all
	if(len(sequenceIdx) == 0):
		return

	# Check that sequenceIdx, sequenceLength, replicationRound, and
	# chromosomeIndex are equal length, numOric should be half the
	# size (4 DNAP/fork, only 2 oriC/fork)
	assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(chromosomeIndex) == 4*(numOric - 1))

	## Update polymerases mass to account for already completed DNA
	# Determine the sequences of already-replicated DNA
	sequences = kb.process.replication.replication_sequences
	sequenceElongations = np.array(sequenceLength, dtype=np.int64)
	massIncreaseDna = computeMassIncrease(
			np.tile(sequences,(len(sequenceIdx) / 4,1)),
			sequenceElongations,
			kb.process.replication.replicationMonomerWeights.asNumber(units.fg)	
			)

	# Update the attributes of replicating DNA polymerases
	oricCenter = kb.constants
	dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', len(sequenceIdx))
	dnaPoly.attrIs(
		sequenceIdx = np.array(sequenceIdx),
		sequenceLength = np.array(sequenceLength),
		replicationRound = np.array(replicationRound),
		chromosomeIndex = np.array(chromosomeIndex),
		massDiff_DNA = massIncreaseDna,
		)

	# Update the attributes of the partially-replicated origins of replication
	oriC = uniqueMolCntr.objectsNew('originOfReplication', numOric)

def setDaughterInitialConditions(sim, kb):
	assert sim._inheritedStatePath != None

	bulk_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "BulkMolecules"))
	sim.states["BulkMolecules"].tableLoad(bulk_table_reader, 0)

	# bulk_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "BulkChromosome"))
	# sim.states["BulkChromosome"].tableLoad(bulk_table_reader, 0)

	unique_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "UniqueMolecules"))
	sim.states["UniqueMolecules"].tableLoad(unique_table_reader, 0)

	time_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "Time"))
	initialTime = TableReader(os.path.join(sim._inheritedStatePath, "Time")).readAttribute("initialTime")
	sim._initialTime = initialTime


def determineChromosomeState(C, D, tau, replication_length):
	"""
	determineChromosomeState

	Purpose: calculates the number and position of replicating DNA polymerases
			 at the beginning of the cell cycle.

	Inputs: C  - the C period of the cell, the length of time between
			replication initiation and replication completion.
			D  - the D period of the cell, the length of time between
			completing replication of the chromosome and division of the cell.
			tau - the doubling time of the cell
			replication_length - the amount of DNA to be replicated per fork,
			usually half of the genome, in base-pairs

	Outputs: a tuple of vectors for input into the dnaPoly.attrIs() function
			of the initializeReplication() function in initial_conditions.py.
			These  vectors are, in order:
			sequenceIdx - an index for each of the four types/directions of DNA
						replication - leading and lagging strand of the forward
						and the reverse fork = 4 total. This vector is always
						simply [0,1,2,3] repeated once for each replication
						event. Ie for three active replication events (6 forks,
						12 polymerases) sequenceIdx = [0,1,2,3,0,1,2,3,0,1,2,3]
			sequenceLength - the position in the genome that each polymerase 
						referenced in sequenceIdx has reached, in base-pairs.
						This is handled such that even though in reality some 
						polymerases are replicating in different directions all
						inputs here are as though each starts at 0 and goes up
						the the number of base-pairs to be replicated.
			replicationRound - a unique integer stating in which replication 
						generation the polymerase referenced by sequenceIdx.
						Each time all origins of replication in the cell fire,
						a new replication generation has started. This array 
						is integer-valued, and counts from 0 (the oldest
						generation) up to n (the most recent initiation) event.
			chromosomeIndex - indicator variable for which daughter cell 
						should inherit which polymerase at division. This
						array is only relevant to draw distinctions within a 
						generaation of replicationRound, now between rounds.
						Within each generation in replicationRound (run of
						numbers with the same value), half should have
						chromosomeIndex = 0, and half chromosomeIndex
						= 1. This should be contiguous halves, ie
						[0,0,0,0,1,1,1,1], NOT interspersed like
						[0,1,0,1,0,1,0,1]. The half-and-half rule is excepted
						for any replication generation with only 4 polymerases
						(the oldest replication generation should be the only
						one	with fewer than 8 polymerases). In this case 
						chromosomeIndex doesn't matter/is effectively NaN,
						but is set to all 0's to prevent conceptually dividing
						a single chromosome between two daughter cells.

	Notes: if NO polymerases are active at the start of the cell cycle,
			equivalent to the C + D periods being shorter than the doubling
			rate tau, then this function returns empty lists. dnaPoly.attrIs()
			should not be run in this case, as no DNA replication will be
			underway.
	"""

	## Error check inputs

	# Check that all inputs have units
	assert (units.hasUnit(C)), 'C must have units'
	assert (units.hasUnit(D)), 'D must have units'
	assert (units.hasUnit(tau)), 'tau must have units'
	assert (units.hasUnit(replication_length)) and (units.getUnit(replication_length).strUnit() == units.nt.strUnit()), 'replication_length must have units of units.nt.'

	# All inputs must be positive numbers
	assert (C.asNumber(units.min) >= 0), "C value can't be negative."
	assert (D.asNumber(units.min) >= 0), "D value can't be negative."
	assert (tau.asNumber(units.min) >= 0), "tau value can't be negative."
	assert (replication_length.asNumber(units.nt) >= 0), "replication_length value can't be negative."

	# Require that D is shorter than tau - time between completing replication
	# and dividing must be shorter than the time between two divisions.
	assert (D.asNumber(units.min) < tau.asNumber(units.min)), 'The D period must be shorter than the doubling time tau.'


	# Number active replication generations (can be many initiations per gen.)
	limit = np.floor((C.asNumber(units.min) + D.asNumber(units.min))/tau.asNumber(units.min))

	# Initialize arrays to be returned
	sequenceIdx = []
	sequenceLength = []
	replicationRound = []
	chromosomeIndex = []

	# Loop through the generations of replication which are active (if limit = 0
	# skips loop entirely --> no active replication generations)
	n = 1;
	while n <= limit:
		# Determine at what base each strand of a given replication event should start
		# Replication forks should be at base (1 - (n*tau - D)/(C))(basepairs to be replicated)

		ratio = (1 - ((n*tau - D)/(C)))
		ratio = units.convertNoUnitToNumber(ratio)
		fork_location = np.floor(ratio*(replication_length.asNumber(units.nt)))

		# Add 2^(n-1) replication events (two forks, four strands per inintiaion event)
		num_events = 2 ** (n-1)
		# sequenceIdx refers to the type of elongation - ie forward and
		# reverse, lagging and leading strands.
		sequenceIdx += [0,1,2,3]*num_events
		# sequenceLength refers to how far along in the replication process 
		# this event already is - at what basepair currently. All four are
		# assumed to be equally far along.
		sequenceLength += [fork_location]*4*num_events
		# replicationRound is defined the same way as n - each time a
		# replication round starts, all origins in the cell fire, and are then
		# of the same generation and round number
		replicationRound += [(n-1)]*4*num_events
		# WITHIN each round, chromosomeIndex uniquely identifies individual
		# origin initaion points. Loop through each intiation event in this 
		# generation (2 forks, 4 polymerases each), assign it an increaing,
		# unique number, starting at zero.
		chromosomeIndex += [0]*2*num_events + [1]*2*num_events

		n += 1

	# The first replication generation should not be divided, so set all values
	# to 0 (effectively NaN, the first four values are not used)
	if len(chromosomeIndex):
		chromosomeIndex[:4] = [0,0,0,0]

	return (sequenceIdx, sequenceLength, replicationRound, chromosomeIndex)



def determineNumOriC(C, D, tau):
	"""
	determineNumOriC

	Purpose: calculates the number of OriC's in a cell upon initiation, determined by the replication state of the chromosome.

	Inputs: C  - the C period of the cell, the length of time between
			replication initiation and replication completion.
			D  - the D period of the cell, the length of time between
			completing replication of the chromosome and division of the cell.
			tau - the doubling time of the cell
			replication_length - the amount of DNA to be replicated per fork,
			usually half of the genome, in base-pairs

	Outputs: the number of OriC's in the cell at initiation.	
	"""

	# Number active replication generations (can be many initiations per gen.)
	total_active_initiations = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())

	return int(2 ** (total_active_initiations))
