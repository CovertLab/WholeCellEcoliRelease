"""
divide_cell.py
Functions needed for division from mother cell into daughter cells
"""
from __future__ import division

import numpy as np
import os
from copy import deepcopy

from wholecell.io.tablewriter import TableWriter

from wholecell.utils import units

BINOMIAL_COEFF = 0.5

def divide_cell(sim):
	# Assign data from simulation required
	randomState = sim.randomState

	bulkMolecules = sim.states['BulkMolecules']
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
	d1_uniqueMolCntr, d2_uniqueMolCntr, daughter_elng_rates = divideUniqueMolecules(uniqueMolecules, randomState, chromosome_counts, sim)

	# Save divded containers
	saveContainer(d1_bulkMolCntr, os.path.join(sim._outputDir, "Daughter1", "BulkMolecules"))
	saveContainer(d2_bulkMolCntr, os.path.join(sim._outputDir, "Daughter2", "BulkMolecules"))
	saveContainer(d1_uniqueMolCntr, os.path.join(sim._outputDir, "Daughter1", "UniqueMolecules"))
	saveContainer(d2_uniqueMolCntr, os.path.join(sim._outputDir, "Daughter2", "UniqueMolecules"))

	import cPickle
	cPickle.dump(daughter_elng_rates["d1_elng_rate"], open(os.path.join(sim._outputDir, "Daughter1", "ElngRate.cPickle"),'wb'))
	cPickle.dump(daughter_elng_rates["d2_elng_rate"], open(os.path.join(sim._outputDir, "Daughter2", "ElngRate.cPickle"),'wb'))

	cPickle.dump(daughter_elng_rates["d1_elng_rate_factor"], open(os.path.join(sim._outputDir, "Daughter1", "elng_rate_factor.cPickle"),'wb'))
	cPickle.dump(daughter_elng_rates["d2_elng_rate_factor"], open(os.path.join(sim._outputDir, "Daughter2", "elng_rate_factor.cPickle"),'wb'))



	# Save daughter cell initial time steps
	saveTime(sim.time(), os.path.join(sim._outputDir, "Daughter1", "Time"), sim.timeStepSec())
	saveTime(sim.time(), os.path.join(sim._outputDir, "Daughter2", "Time"), sim.timeStepSec())

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

	# Divide evenly between both daughters if even number of chromosomes or give an extra to one cell if odd number
	if full_chromosome_count % 2 == 0:
		d1_chromosome_count = full_chromosome_count // 2
		d2_chromosome_count = full_chromosome_count // 2
	else:
		d1_chromosome_count = full_chromosome_count // 2
		d1_chromosome_count += randomState.binomial(1, p = BINOMIAL_COEFF)
		d2_chromosome_count = full_chromosome_count - d1_chromosome_count

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

	## Set basal transcription values to 1
	# TODO: Add assertion that these elements have zero mass
	d1_bulk_molecules_container.countsIs(1, bulkMolecules.divisionIds['setTo1'])
	d2_bulk_molecules_container.countsIs(1, bulkMolecules.divisionIds['setTo1'])

	return d1_bulk_molecules_container, d2_bulk_molecules_container

def divideUniqueMolecules(uniqueMolecules, randomState, chromosome_counts, sim):
	d1_unique_molecules_container = uniqueMolecules.container.emptyLike()
	d2_unique_molecules_container = uniqueMolecules.container.emptyLike()

	uniqueMoleculesToDivide = deepcopy(uniqueMolecules.uniqueMoleculeDefinitions)

	# Divide unique molecules binomially
	for moleculeName, moleculeAttributeDict in uniqueMoleculesToDivide.iteritems():
		if moleculeName == 'dnaPolymerase' or moleculeName == 'originOfReplication' or moleculeName == 'fullChromosome' or moleculeName == 'activeRibosome':# or moleculeName == 'activeRnaPoly':# or moleculeName == 'activeRibosome':
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

	# Unequally divide active RNAP
	# moleculeSet = uniqueMolecules.container.objectsInCollection('activeRnaPoly')
	# moleculeAttributeDict = uniqueMoleculesToDivide['activeRnaPoly']
	# if len(moleculeSet) > 0:

	# 	# coeff = randomState.normal(1.0, 0.6) * 0.5
	# 	# n_d1 = randomState.binomial(len(moleculeSet), p = coeff)

	# 	# n_d1 = int(np.round(len(moleculeSet) * randomState.normal(1.0, 0.3) * 0.5))
	# 	# n_d2 = len(moleculeSet) - n_d1

	# 	mean = len(moleculeSet) * 0.5
	# 	sd = 0.2 * mean
	# 	n_d1 = int(np.round(np.random.normal(mean, sd)))

	# 	n_d2 = len(moleculeSet) - n_d1
	# 	assert n_d1 + n_d2 == len(moleculeSet)

	# 	d1_bool = np.zeros(len(moleculeSet), dtype = bool)
	# 	d2_bool = np.zeros(len(moleculeSet), dtype = bool)

	# 	d1_indexes = randomState.choice(range(len(moleculeSet)), size = n_d1, replace = False)
	# 	d1_bool[d1_indexes] = True
	# 	d2_bool = np.logical_not(d1_bool)

	# 	d1_dividedAttributesDict = {}
	# 	d2_dividedAttributesDict = {}
	# 	for moleculeAttribute in moleculeAttributeDict.iterkeys():
	# 		d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
	# 		d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

	# 	d1_unique_molecules_container.objectsNew('activeRnaPoly', n_d1, **d1_dividedAttributesDict)
	# 	d2_unique_molecules_container.objectsNew('activeRnaPoly', n_d2, **d2_dividedAttributesDict)


	# Unequally divide active ribosomes
	moleculeSet = uniqueMolecules.container.objectsInCollection('activeRibosome')
	moleculeAttributeDict = uniqueMoleculesToDivide['activeRibosome']
	if len(moleculeSet) > 0:

		polyElng = sim.processes["PolypeptideElongation"]

		environmentalElongationRate = polyElng.ribosomeElongationRateDict[polyElng.currentNutrients].asNumber(units.aa / units.s)

		elngRate = np.min([polyElng.ribosomeElongationRateDict[polyElng.currentNutrients].asNumber(units.aa / units.s), 21.])
		nRibosomes = len(uniqueMolecules.container.objectsInCollection("activeRibosome"))
		noiseMultiplier = 1.
		if sim._growthRateNoise:
			noiseMultiplier = randomState.normal(1, 0.2)
		translationCapacity = elngRate * nRibosomes * noiseMultiplier

		# mean = len(moleculeSet) * 0.5
		# sd = 0.000001 * mean

		# n_d1 = int(np.round(randomState.normal(mean, sd)))
		n_d1 = randomState.binomial(len(moleculeSet), 0.5)
		n_d2 = len(moleculeSet) - n_d1
		assert n_d1 + n_d2 == len(moleculeSet)

		d1_rib_elng_rate = np.min([(translationCapacity / 2) / n_d1, 21.])
		d2_rib_elng_rate = np.min([(translationCapacity / 2) / n_d2, 21.])

		daughter_elng_rates = {
						"d1_elng_rate" : d1_rib_elng_rate,
						"d2_elng_rate" : d2_rib_elng_rate,
						"d1_elng_rate_factor" : d1_rib_elng_rate / environmentalElongationRate,
						"d2_elng_rate_factor" : d2_rib_elng_rate / environmentalElongationRate,
						}

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

		d1_unique_molecules_container.objectsNew('activeRibosome', n_d1, **d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('activeRibosome', n_d2, **d2_dividedAttributesDict)



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

		print "State at division:"
		sequenceIdx, sequenceLengths, replicationRound, chromosomeIndex = moleculeSet.attrs(
			'sequenceIdx', 'sequenceLength', 'replicationRound', 'chromosomeIndex'
			)
		print "sequenceIdx: {}".format(sequenceIdx)
		print "sequenceLengths: {}".format(sequenceLengths)
		print "replicationRound: {}".format(replicationRound)
		print "chromosomeIndex: {}".format(chromosomeIndex)

		if d1_chromosome_count + d2_chromosome_count < 2:
			d1_bool = np.zeros(len(moleculeSet), dtype = bool)
			d2_bool = np.zeros(len(moleculeSet), dtype = bool)

			d1_bool[:] = d1_chromosome_count
			d2_bool[:] = d2_chromosome_count
		else:
			d1_bool = np.zeros(len(moleculeSet), dtype = bool)
			d2_bool = np.zeros(len(moleculeSet), dtype = bool)

			# give polymerases associated with extra chromosome to daughter that gets the extra chromosome
			if d1_chromosome_count < d2_chromosome_count:
				d1Index = 1
				d2Index = 0
			else:
				d1Index = 0
				d2Index = 1

			for roundIdx in np.unique(replicationRound):
				# replicationRound indexes all dna polymerases started at the same time across all oriC
				# chromosomeIndex indexes all dna polymerases started at the same time at EACH oriC (i.e. one should go to d1, one to d2)
				d1_bool = np.logical_or(
					np.logical_and(replicationRound == roundIdx, chromosomeIndex % 2 == d1Index),
					d1_bool
					)
				d2_bool = np.logical_or(
					np.logical_and(replicationRound == roundIdx, chromosomeIndex % 2 == d2Index),
					d2_bool
					)

		print "d1_bool: {}".format(d1_bool)
		print "d2_bool: {}".format(d2_bool)

		n_d1 = d1_bool.sum()
		n_d2 = d2_bool.sum()

		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
			d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

		for replicationRound in np.unique(d2_dividedAttributesDict['replicationRound']):
			replicationRoundIndexes = d2_dividedAttributesDict['replicationRound'] == replicationRound
			replicationForks = int(np.ceil(replicationRoundIndexes.sum() / 4))
			if replicationForks >= 2:
				for index in [0,1,2,3]:
					# if not possible to have uneven number of partial chromosomes
					num_index = d2_dividedAttributesDict['sequenceIdx'][replicationRoundIndexes] == index
					new_value = np.zeros(num_index.sum())
					for fork in range(replicationForks):
						new_value[fork * new_value.size / replicationForks:(fork + 1) * new_value.size / replicationForks] = fork

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

		d1_replicationForks = int(np.ceil(len(d1_unique_molecules_container.objectsInCollection('dnaPolymerase')) / 4))
		n_d1 = d1_chromosome_count + d1_replicationForks
		d1_bool[:n_d1] = 1
		d2_bool = np.logical_not(d1_bool)
		n_d2 = d2_bool.sum()

		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
			d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

		d1_unique_molecules_container.objectsNew('originOfReplication', n_d1, **d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('originOfReplication', n_d2, **d2_dividedAttributesDict)

	# Divide fullChromosomes with chromosome
	moleculeSet = uniqueMolecules.container.objectsInCollection('fullChromosome')
	moleculeAttributeDict = uniqueMoleculesToDivide['fullChromosome']
	if len(moleculeSet) > 0:
		d1_chromosome_count = chromosome_counts['d1_chromosome_count']
		d2_chromosome_count = chromosome_counts['d2_chromosome_count']

		d1_bool = np.ones(len(moleculeSet), dtype = bool)
		d2_bool = np.ones(len(moleculeSet), dtype = bool)

		n_d1 = len(moleculeSet)
		n_d2 = len(moleculeSet)

		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
			d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

		d1_unique_molecules_container.objectsNew('fullChromosome', n_d1, **d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('fullChromosome', n_d2, **d2_dividedAttributesDict)

		print "grep_marker divide cell - fullChromosome to daughter 1: {}".format(n_d1)
		print "grep_marker divide cell - fullChromosome to daughter 2: {}".format(n_d2)

	return d1_unique_molecules_container, d2_unique_molecules_container, daughter_elng_rates

def saveContainer(container, path):
	table_writer = TableWriter(path)
	container.tableCreate(table_writer)
	container.tableAppend(table_writer)

def saveTime(finalTime, path, timeStepSec):
	timeFile = TableWriter(path)

	# Metadata
	timeFile.writeAttributes(
		initialTime = finalTime + timeStepSec
		)

	timeFile.close()
