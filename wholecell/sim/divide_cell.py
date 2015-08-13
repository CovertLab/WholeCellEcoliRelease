"""
divide_cell.py
Functions needed for division from mother cell into daughter cells
"""
from __future__ import division

import numpy as np
import os
from copy import deepcopy

from wholecell.io.tablewriter import TableWriter

BINOMIAL_COEFF = 0.5

def divide_cell(sim):
	# Assign data from simulation required
	randomState = sim.randomState

	bulkMolecules = sim.states['BulkMolecules']
	bulkChromosome = sim.states['BulkChromosome']
	uniqueMolecules = sim.states['UniqueMolecules']

	# Create output directories
	try:
		os.mkdir(os.path.join(sim._outputDir, "Daughter1"))
	except OSError:
		pass

	try:
		os.mkdir(os.path.join(sim._outputDir, "Daughter2"))
	except OSError:
		pass

	# Calculate chromosome division
	# Used in dividing both bulk molecules and unique molecules
	chromosome_counts = chromosomeDivision(bulkMolecules, randomState)

	# Create divded containers
	d1_bulkMolCntr, d2_bulkMolCntr = divideBulkMolecules(bulkMolecules, randomState, chromosome_counts)
	# d1_bulkChrmCntr, d2_bulkChrmCntr = divideBulkChromosome(bulkChromosome, randomState)
	d1_uniqueMolCntr, d2_uniqueMolCntr = divideUniqueMolecules(uniqueMolecules, randomState, chromosome_counts)

	# Save divded containers
	saveContainer(d1_bulkMolCntr, os.path.join(sim._outputDir, "Daughter1", "BulkMolecules"))
	saveContainer(d2_bulkMolCntr, os.path.join(sim._outputDir, "Daughter2", "BulkMolecules"))
	# saveContainer(d1_bulkChrmCntr, os.path.join(sim._outputDir, "Daughter1", "BulkChromosome"))
	# saveContainer(d2_bulkChrmCntr, os.path.join(sim._outputDir, "Daughter2", "BulkChromosome"))
	saveContainer(d1_uniqueMolCntr, os.path.join(sim._outputDir, "Daughter1", "UniqueMolecules"))
	saveContainer(d2_uniqueMolCntr, os.path.join(sim._outputDir, "Daughter2", "UniqueMolecules"))

	# Save daughter cell initial time steps
	saveTime(sim.time(), os.path.join(sim._outputDir, "Daughter1", "Time"))
	saveTime(sim.time(), os.path.join(sim._outputDir, "Daughter2", "Time"))

def chromosomeDivision(bulkMolecules, randomState):
	partial_chromosome_counts = bulkMolecules.container.counts(bulkMolecules.divisionIds['partialChromosome'])

	uneven_counts = partial_chromosome_counts - partial_chromosome_counts.min()
	if uneven_counts.any():
		raise Exception("You won the lottery! There is an uneven number of partial chromosomes...")

	# Transform any leftover partial chromosomes into full chromosome for convienence
	# this should have happened in simulation process but in this simulation we got
	# lucky and missed this before the final time-step.
	bulkMolecules.container.countInc(
		partial_chromosome_counts.min(),
		bulkMolecules.divisionIds['fullChromosome'][0]
		)
	full_chromosome_count = bulkMolecules.container.count(bulkMolecules.divisionIds['fullChromosome'][0])

	if full_chromosome_count == 1:
		d1_chromosome_count = randomState.binomial(full_chromosome_count, p = BINOMIAL_COEFF)
		d2_chromosome_count = full_chromosome_count - d1_chromosome_count
	elif full_chromosome_count == 2:
		d1_chromosome_count = 1
		d2_chromosome_count = 1
	else:
		raise Exception("Un-accounted for number of chromosomes at division!")

	return {"d1_chromosome_count" : d1_chromosome_count, "d2_chromosome_count" : d2_chromosome_count}

def divideBulkMolecules(bulkMolecules, randomState, chromosome_counts):
	d1_bulk_molecules_container = bulkMolecules.container.emptyLike()
	d2_bulk_molecules_container = bulkMolecules.container.emptyLike()

	## Divide most bulk molecules binomially
	molecule_counts = bulkMolecules.container.counts(bulkMolecules.divisionIds['binomial'])

	d1_counts = randomState.binomial(molecule_counts, p = BINOMIAL_COEFF)
	d2_counts = molecule_counts - d1_counts

	assert np.all(d1_counts + d2_counts == molecule_counts)

	d1_bulk_molecules_container.countsIs(d1_counts, bulkMolecules.divisionIds['binomial'])
	d2_bulk_molecules_container.countsIs(d2_counts, bulkMolecules.divisionIds['binomial'])

	## NOTE: THERE ARE NO EQUAL DIVISION MOLECULES RIGHT NOW SO THIS IS COMMENTED
	## Handle molecules that should be divded equally
	## If there are an odd number of molecules they are stochastically allocated.
	# molecule_counts = bulkMolecules.container.counts(bulkMolecules.divisionIds['equally'])

	# uneven_counts = molecule_counts % 2
	# d1_uneven_counts = randomState.binomial(uneven_counts, p = BINOMIAL_COEFF)
	# d2_uneven_counts = uneven_counts - d1_uneven_counts

	# d1_counts = molecule_counts / 2 + d1_uneven_counts
	# d2_counts = molecule_counts / 2 + d2_uneven_counts

	# assert all(d1_counts + d2_counts == molecule_counts)

	# d1_bulk_molecules_container.countsIs(d1_counts, bulkMolecules.divisionIds['equally'])
	# d2_bulk_molecules_container.countsIs(d2_counts, bulkMolecules.divisionIds['equally'])

	## Handle chromosome division
	full_chromosome_count = bulkMolecules.container.count(bulkMolecules.divisionIds['fullChromosome'][0])

	d1_chromosome_count = chromosome_counts['d1_chromosome_count']
	d2_chromosome_count = chromosome_counts['d2_chromosome_count']

	assert np.all(d1_chromosome_count + d2_chromosome_count == full_chromosome_count)

	d1_bulk_molecules_container.countIs(d1_chromosome_count, bulkMolecules.divisionIds['fullChromosome'][0])
	d2_bulk_molecules_container.countIs(d2_chromosome_count, bulkMolecules.divisionIds['fullChromosome'][0])

	return d1_bulk_molecules_container, d2_bulk_molecules_container

def divideBulkChromosome(bulkChromosome, randomState):
	chromosome_location_counts = bulkChromosome.container.counts()
	d1_bulk_chromosome_container = bulkChromosome.container.emptyLike()
	d2_bulk_chromosome_container = bulkChromosome.container.emptyLike()

	# TODO: Delete bulkChromosome entirely!
	# if dnaReplicationComplete:
	# 	assert all(chromosome_location_counts == 2)

	# 	d1_counts = chromosome_location_counts / 2
	# 	d2_counts = chromosome_location_counts / 2

	# 	assert all(d1_counts + d2_counts == chromosome_location_counts)

	# 	d1_bulk_chromosome_container.countsIs(d1_counts)
	# 	d2_bulk_chromosome_container.countsIs(d2_counts)
	# elif chromosomeToDaughter1:
	# 	d1_bulk_chromosome_container.countsIs(chromosome_location_counts)
	# else:
	# 	d2_bulk_chromosome_container.countsIs(chromosome_location_counts)

	return d1_bulk_chromosome_container, d2_bulk_chromosome_container

def divideUniqueMolecules(uniqueMolecules, randomState, chromosome_counts):
	d1_unique_molecules_container = uniqueMolecules.container.emptyLike()
	d2_unique_molecules_container = uniqueMolecules.container.emptyLike()

	uniqueMoleculesToDivide = deepcopy(uniqueMolecules.uniqueMoleculeDefinitions)

	# Divide unique molecules binomially
	for moleculeName, moleculeAttributeDict in uniqueMoleculesToDivide.iteritems():
		if moleculeName == 'dnaPolymerase' or moleculeName == 'originOfReplication':
			# NOTE: We are not dividing dna polymerase binomially!
			continue

		# Get set of molecules to divide and calculate number going to daugher one and daughter two
		moleculeSet = uniqueMolecules.container.objectsInCollection(moleculeName)
		if len(moleculeSet) > 0:
			n_d1 = randomState.binomial(len(moleculeSet), p = BINOMIAL_COEFF)
			n_d2 = len(moleculeSet) - n_d1
			assert n_d1 + n_d2 == len(moleculeSet)

			# Randomly boolean index molecules in mother such that each daugher gets amount calculated above
			d1_bool = np.zeros(len(moleculeSet), dtype = bool)
			d2_bool = np.zeros(len(moleculeSet), dtype = bool)
			d1_indexes = randomState.choice(range(len(moleculeSet)), size = n_d1, replace = False)
			d1_bool[d1_indexes] = True
			d2_bool = np.logical_not(d1_bool)

			d1_dividedAttributesDict = {}
			d2_dividedAttributesDict = {}
			for moleculeAttribute in moleculeAttributeDict.iterkeys():
				d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
				d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

			d1_unique_molecules_container.objectsNew(moleculeName, n_d1, **d1_dividedAttributesDict)
			d2_unique_molecules_container.objectsNew(moleculeName, n_d2, **d2_dividedAttributesDict)

	# Divide dna polymerase with chromosome
	# Get set of molecules to divide and calculate number going to daugher one and daughter two
	moleculeSet = uniqueMolecules.container.objectsInCollection('dnaPolymerase')
	moleculeAttributeDict = uniqueMoleculesToDivide['dnaPolymerase']
	if len(moleculeSet) > 0:
		d1_chromosome_count = chromosome_counts['d1_chromosome_count']
		d2_chromosome_count = chromosome_counts['d2_chromosome_count']

		sequenceIdx, sequenceLengths, replicationRound, chromosomeIndex = moleculeSet.attrs(
			'sequenceIdx', 'sequenceLength', 'replicationRound', 'chromosomeIndex'
			)

		if d1_chromosome_count + d2_chromosome_count < 2:
			d1_bool = np.zeros(len(moleculeSet), dtype = bool)
			d2_bool = np.zeros(len(moleculeSet), dtype = bool)

			d1_bool[:] = d1_chromosome_count
			d2_bool[:] = d2_chromosome_count
		elif d1_chromosome_count + d2_chromosome_count == 2:
			d1_bool = np.zeros(len(moleculeSet), dtype = bool)
			d2_bool = np.zeros(len(moleculeSet), dtype = bool)
			for roundIdx in np.unique(replicationRound):
				# replicationRound indexes all dna polymerases started at the same time across all oriC
				# chromosomeIndex indexes all dna polymerases started at the same time at EACH oriC (i.e. one should go to d1, one to d2)
				d1_bool = np.logical_or(
					np.logical_and(replicationRound == roundIdx, chromosomeIndex == 0),
					d1_bool
					)
				d2_bool = np.logical_or(
					np.logical_and(replicationRound == roundIdx, chromosomeIndex == 1),
					d2_bool
					)
		else:
			raise Exception("Too may chromosomes!")
				
		n_d1 = d1_bool.sum()
		n_d2 = d2_bool.sum()

		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
			d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

		# Reset chromosomeIndex for next cell cycle in all forks that are divided
		for replicationRound in np.unique(d1_dividedAttributesDict['replicationRound']):
			replicationRoundIndexes = d1_dividedAttributesDict['replicationRound'] == replicationRound
			if replicationRoundIndexes.sum() >= 8:
				for index in [0,1,2,3]:
					num_index = d1_dividedAttributesDict['sequenceIdx'][replicationRoundIndexes] == index
					new_value = np.zeros(num_index.sum())
					new_value[:new_value.size / 2] = 1
					d1_dividedAttributesDict['chromosomeIndex'][num_index] = new_value

				# zero_index = d1_dividedAttributesDict['sequenceIdx'][replicationRoundIndexes] == 0
				# one_index = d1_dividedAttributesDict['sequenceIdx'][replicationRoundIndexes] == 1
				# two_index = d1_dividedAttributesDict['sequenceIdx'][replicationRoundIndexes] == 2
				# three_index = d1_dividedAttributesDict['sequenceIdx'][replicationRoundIndexes] == 3

				# new_value = np.zeros(zero_index.sum())
				# new_value[:new_value.size / 2] = 1
				# d1_dividedAttributesDict['chromosomeIndex'][zero_index] = new_value

		for replicationRound in np.unique(d2_dividedAttributesDict['replicationRound']):
			replicationRoundIndexes = d2_dividedAttributesDict['replicationRound'] == replicationRound
			if replicationRoundIndexes.sum() >= 8:
				for index in [0,1,2,3]:
					num_index = d2_dividedAttributesDict['sequenceIdx'][replicationRoundIndexes] == index
					new_value = np.zeros(num_index.sum())
					new_value[:new_value.size / 2] = 0
					d2_dividedAttributesDict['chromosomeIndex'][num_index] = new_value

		d1_unique_molecules_container.objectsNew('dnaPolymerase', n_d1, **d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('dnaPolymerase', n_d2, **d2_dividedAttributesDict)

		# Divide oriCs with chromosome

	moleculeSet = uniqueMolecules.container.objectsInCollection('originOfReplication')
	moleculeAttributeDict = uniqueMoleculesToDivide['originOfReplication']
	if len(moleculeSet) > 0:
		d1_chromosome_count = chromosome_counts['d1_chromosome_count']
		d2_chromosome_count = chromosome_counts['d2_chromosome_count']

		d1_bool = np.zeros(len(moleculeSet), dtype = bool)
		d2_bool = np.zeros(len(moleculeSet), dtype = bool)
		if d1_chromosome_count + d2_chromosome_count < 2:
			d1_bool[:] = d1_chromosome_count
			d2_bool[:] = d2_chromosome_count

		elif d1_chromosome_count + d2_chromosome_count == 2:
			d1_bool[:len(moleculeSet) / 2.] = 1
			d2_bool = np.logical_not(d1_bool)

		else:
			raise Exception("Too may chromosomes!")
				
		n_d1 = d1_bool.sum()
		n_d2 = d2_bool.sum()

		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
			d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

		d1_unique_molecules_container.objectsNew('originOfReplication', n_d1, **d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('originOfReplication', n_d2, **d2_dividedAttributesDict)


	return d1_unique_molecules_container, d2_unique_molecules_container

def saveContainer(container, path):
	table_writer = TableWriter(path)
	container.tableCreate(table_writer)
	container.tableAppend(table_writer)

def saveTime(finalTime, path):
	timeFile = TableWriter(path)

	# Metadata
	timeFile.writeAttributes(
		initialTime = finalTime + 1
		)

	timeFile.close()
