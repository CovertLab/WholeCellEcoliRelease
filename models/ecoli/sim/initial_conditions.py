
"""

TODO:
- document math
- replace fake metabolite concentration targets with measured metabolite concentration targets
- raise/warn if physiological metabolite concentration targets appear to be smaller than what
 is needed at this time step size

"""

from __future__ import division

from itertools import izip

import numpy as np
import os

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.utils.fitting import normalize, countsFromMassAndExpression, calcProteinCounts, massesAndCountsToAddForHomeostaticTargets
from wholecell.utils.polymerize import buildSequences, computeMassIncrease
from wholecell.utils import units



from wholecell.io.tablereader import TableReader

def calcInitialConditions(sim, sim_data):
	assert sim._inheritedStatePath == None
	randomState = sim.randomState

	massCoeff = 1.0
	if sim._massDistribution:
		massCoeff = randomState.normal(loc = 1.0, scale = 0.1)

	bulkMolCntr = sim.states['BulkMolecules'].container
	uniqueMolCntr = sim.states["UniqueMolecules"].container

	# Set up states
	initializeBulkMolecules(bulkMolCntr, sim_data, randomState, massCoeff)
	initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

def initializeBulkMolecules(bulkMolCntr, sim_data, randomState, massCoeff):

	## Set protein counts from expression
	initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff)

	## Set RNA counts from expression
	initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff)

	## Set DNA
	initializeDNA(bulkMolCntr, sim_data, randomState)

	## Set other biomass components
	initializeSmallMolecules(bulkMolCntr, sim_data, randomState, massCoeff)

	## Set constitutive expression
	initializeConstitutiveExpression(bulkMolCntr, sim_data, randomState)

def initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState):
	initializeReplication(bulkMolCntr, uniqueMolCntr, sim_data)

def initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff):

	monomersView = bulkMolCntr.countsView(sim_data.process.translation.monomerData["id"])
	monomerMass = massCoeff * sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])["proteinMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	# TODO: unify this logic with the fitter so it doesn't fall out of step
	# again (look at the calcProteinCounts function)

	monomerExpression = normalize(
		sim_data.process.transcription.rnaExpression[sim_data.condition][sim_data.relation.rnaIndexToMonomerMapping] *
		sim_data.process.translation.translationEfficienciesByMonomer /
		(np.log(2) / sim_data.conditionToDoublingTime[sim_data.condition].asNumber(units.s) + sim_data.process.translation.monomerData["degRate"].asNumber(1 / units.s))
		)

	nMonomers = countsFromMassAndExpression(
		monomerMass.asNumber(units.g),
		sim_data.process.translation.monomerData["mw"].asNumber(units.g/units.mol),
		monomerExpression,
		sim_data.constants.nAvogadro.asNumber(1/units.mol)
		)

	monomersView.countsIs(
		randomState.multinomial(nMonomers, monomerExpression)
		)

def initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff):

	rnaView = bulkMolCntr.countsView(sim_data.process.transcription.rnaData["id"])
	rnaMass = massCoeff * sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor

	rnaExpression = normalize(sim_data.process.transcription.rnaExpression[sim_data.condition])

	nRnas = countsFromMassAndExpression(
		rnaMass.asNumber(units.g),
		sim_data.process.transcription.rnaData["mw"].asNumber(units.g/units.mol),
		rnaExpression,
		sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		)

	rnaView.countsIs(
		randomState.multinomial(nRnas, rnaExpression)
		)

def initializeDNA(bulkMolCntr, sim_data, randomState):

	chromosomeView = bulkMolCntr.countsView(sim_data.moleculeGroups.fullChromosome)
	chromosomeView.countsIs([1])

# TODO: remove checks for zero concentrations (change to assertion)
# TODO: move any rescaling logic to KB/fitting
def initializeSmallMolecules(bulkMolCntr, sim_data, randomState, massCoeff):
	avgCellFractionMass = sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])

	mass = massCoeff * (avgCellFractionMass["proteinMass"] + avgCellFractionMass["rnaMass"] + avgCellFractionMass["dnaMass"]) / sim_data.mass.avgCellToInitialCellConvFactor

	concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
		sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel][0][1]
		)
	concDict.update(sim_data.mass.getBiomassAsConcentrations(sim_data.conditionToDoublingTime[sim_data.condition]))
	moleculeIds = sorted(concDict)
	moleculeConcentrations = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in moleculeIds])

	massesToAdd, countsToAdd = massesAndCountsToAddForHomeostaticTargets(
		mass,
		moleculeIds,
		moleculeConcentrations,
		sim_data.getter.getMass(moleculeIds),
		sim_data.constants.cellDensity,
		sim_data.constants.nAvogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		moleculeIds
		)

def initializeConstitutiveExpression(bulkMolCntr, sim_data, randomState):
	recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
	alphaNames = [x for x in recruitmentColNames if x.endswith("__alpha")]
	alphaView = bulkMolCntr.countsView(alphaNames)
	alphaView.countsIs(1)

def initializeReplication(bulkMolCntr, uniqueMolCntr, sim_data):
	"""
	initializeReplication

	Purpose: Create the appropriate number of replication forks given the cell growth rate.
	"""

	## Determine the number and location of replication forks at the start of the cell cycle
	# Find growth rate constants
	C = sim_data.growthRateParameters.c_period.asUnit(units.min)
	D = sim_data.growthRateParameters.d_period.asUnit(units.min)
	tau = sim_data.conditionToDoublingTime[sim_data.condition].asUnit(units.min)
	genome_length = sim_data.process.replication.genome_length
	replication_length = np.ceil(.5*genome_length) * units.nt

	# Generate arrays specifying appropriate replication conditions
	sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

	# Determine the number of OriC's currently present in the cell
	numOric = determineNumOriC(C, D, tau)
	oriC = uniqueMolCntr.objectsNew('originOfReplication', numOric)

	# Check that sequenceIdx, sequenceLength, replicationRound, and
	# chromosomeIndex are equal length, numOric should be half the
	# size (4 DNAP/fork, only 2 oriC/fork)
	assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(chromosomeIndex) == 4*(numOric - 1))

	# Return if no replication is occuring at all
	if(len(sequenceIdx) == 0):
		return

	## Update polymerases mass to account for already completed DNA
	# Determine the sequences of already-replicated DNA
	sequences = sim_data.process.replication.replication_sequences
	sequenceElongations = np.array(sequenceLength, dtype=np.int64)
	massIncreaseDna = computeMassIncrease(
			np.tile(sequences,(len(sequenceIdx) / 4,1)),
			sequenceElongations,
			sim_data.process.replication.replicationMonomerWeights.asNumber(units.fg)	
			)

	# Update the attributes of replicating DNA polymerases
	oricCenter = sim_data.constants
	dnaPoly = uniqueMolCntr.objectsNew('dnaPolymerase', len(sequenceIdx))
	dnaPoly.attrIs(
		sequenceIdx = np.array(sequenceIdx),
		sequenceLength = np.array(sequenceLength),
		replicationRound = np.array(replicationRound),
		chromosomeIndex = np.array(chromosomeIndex),
		massDiff_DNA = massIncreaseDna,
		)

def setDaughterInitialConditions(sim, sim_data):
	assert sim._inheritedStatePath != None

	bulk_table_reader = TableReader(os.path.join(sim._inheritedStatePath, "BulkMolecules"))
	sim.states["BulkMolecules"].tableLoad(bulk_table_reader, 0)

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
		chromosomeIndex += [0]*2*num_events + [0]*2*num_events

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
