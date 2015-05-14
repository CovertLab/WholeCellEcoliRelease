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
	os.mkdir(os.path.join(sim._outputDir, "Daughter1"))
	os.mkdir(os.path.join(sim._outputDir, "Daughter2"))

	# Determine where chromosome associated components go
	# if chromosome has not segregated
	dnaReplicationComplete = sim._dnaReplicationComplete
	chromosomeToDaughter1 = None
	if not dnaReplicationComplete:
		chromosomeToDaughter1 = False
		if randomState.binomial(1, p = BINOMIAL_COEFF) == 0:
			chromosomeToDaughter1 = True

	# Create divded containers
	d1_bulkMolCntr, d2_bulkMolCntr = divideBulkMolecules(bulkMolecules, randomState, dnaReplicationComplete, chromosomeToDaughter1)
	d1_bulkChrmCntr, d2_bulkChrmCntr = divideBulkChromosome(bulkChromosome, randomState, dnaReplicationComplete, chromosomeToDaughter1)
	d1_uniqueMolCntr, d2_uniqueMolCntr = divideUniqueMolecules(uniqueMolecules, randomState, dnaReplicationComplete, chromosomeToDaughter1)

	# Save divded containers
	saveContainer(d1_bulkMolCntr, os.path.join(sim._outputDir, "Daughter1", "BulkMolecules"))
	saveContainer(d2_bulkMolCntr, os.path.join(sim._outputDir, "Daughter2", "BulkMolecules"))
	saveContainer(d1_bulkChrmCntr, os.path.join(sim._outputDir, "Daughter1", "BulkChromosome"))
	saveContainer(d2_bulkChrmCntr, os.path.join(sim._outputDir, "Daughter2", "BulkChromosome"))
	saveContainer(d1_uniqueMolCntr, os.path.join(sim._outputDir, "Daughter1", "UniqueMolecules"))
	saveContainer(d2_uniqueMolCntr, os.path.join(sim._outputDir, "Daughter2", "UniqueMolecules"))

	# Save daughter cell initial time steps
	saveTime(sim.time(), os.path.join(sim._outputDir, "Daughter1", "Time"))
	saveTime(sim.time(), os.path.join(sim._outputDir, "Daughter2", "Time"))

def divideBulkMolecules(bulkMolecules, randomState, dnaReplicationComplete, chromosomeToDaughter1):
	d1_bulk_molecules_container = bulkMolecules.container.emptyLike()
	d2_bulk_molecules_container = bulkMolecules.container.emptyLike()

	# Divide most bulk molecules binomially

	molecule_counts = bulkMolecules.container.counts(bulkMolecules.divisionIds['binomial'])

	d1_counts = randomState.binomial(molecule_counts, p = BINOMIAL_COEFF)
	d2_counts = molecule_counts - d1_counts

	assert all(d1_counts + d2_counts == molecule_counts)

	d1_bulk_molecules_container.countsIs(d1_counts, bulkMolecules.divisionIds['binomial'])
	d2_bulk_molecules_container.countsIs(d2_counts, bulkMolecules.divisionIds['binomial'])

	# Handle special cases for chromosome separation. If chromosome is decatinated then all
	# chromosome associated bulk molecules will divide equally. Otherwise will parition with
	# one or the other daughter.

	if dnaReplicationComplete:
		molecule_counts = bulkMolecules.container.counts(bulkMolecules.divisionIds['equally'])

		d1_counts = molecule_counts / 2
		d2_counts = molecule_counts / 2

		assert all(d1_counts + d2_counts == molecule_counts)

		d1_bulk_molecules_container.countsIs(d1_counts, bulkMolecules.divisionIds['equally'])
		d2_bulk_molecules_container.countsIs(d2_counts, bulkMolecules.divisionIds['equally'])
	elif chromosomeToDaughter1:
		molecule_counts = bulkMolecules.container.counts(bulkMolecules.divisionIds['with_chromosome'])
		d1_bulk_molecules_container.countsIs(molecule_counts, bulkMolecules.divisionIds['with_chromosome'])
	else:
		molecule_counts = bulkMolecules.container.counts(bulkMolecules.divisionIds['with_chromosome'])
		d2_bulk_molecules_container.countsIs(molecule_counts, bulkMolecules.divisionIds['with_chromosome'])

	return d1_bulk_molecules_container, d2_bulk_molecules_container

def divideBulkChromosome(bulkChromosome, randomState, dnaReplicationComplete, chromosomeToDaughter1):
	chromosome_location_counts = bulkChromosome.container.counts()
	d1_bulk_chromosome_container = bulkChromosome.container.emptyLike()
	d2_bulk_chromosome_container = bulkChromosome.container.emptyLike()

	if dnaReplicationComplete:
		assert all(chromosome_location_counts == 2)

		d1_counts = chromosome_location_counts / 2
		d2_counts = chromosome_location_counts / 2

		assert all(d1_counts + d2_counts == chromosome_location_counts)

		d1_bulk_chromosome_container.countsIs(d1_counts)
		d2_bulk_chromosome_container.countsIs(d2_counts)
	elif chromosomeToDaughter1:
		d1_bulk_chromosome_container.countsIs(chromosome_location_counts)
	else:
		d2_bulk_chromosome_container.countsIs(chromosome_location_counts)

	return d1_bulk_chromosome_container, d2_bulk_chromosome_container

def divideUniqueMolecules(uniqueMolecules, randomState, dnaReplicationComplete, chromosomeToDaughter1):
	d1_unique_molecules_container = uniqueMolecules.container.emptyLike()
	d2_unique_molecules_container = uniqueMolecules.container.emptyLike()

	uniqueMoleculesToDivide = deepcopy(uniqueMolecules.uniqueMoleculeDefinitions)
	if dnaReplicationComplete:
		assert len(uniqueMolecules.container.objectsInCollection('dnaPolymerase')) == 0

	# Divide unique molecules binomially
	for moleculeName, moleculeAttributeDict in uniqueMoleculesToDivide.iteritems():
		if moleculeName == 'dnaPolymerase':
			# NOTE: We are not dividing dna polymerase binomially!
			continue

		# Get set of molecules to divide and calculate number going to daugher one and daughter two
		moleculeSet = uniqueMolecules.container.objectsInCollection(moleculeName)
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
	if not dnaReplicationComplete:
		# Get set of molecules to divide and calculate number going to daugher one and daughter two
		moleculeSet = uniqueMolecules.container.objectsInCollection('dnaPolymerase')
		if len(moleculeSet) == 0:
			raise Exception("No dna polymerase exist even though dna replication has not finished!")

		if chromosomeToDaughter1:
			n_d1 = len(moleculeSet)
		else:
			n_d1 = 0
		n_d2 = len(moleculeSet) - n_d1
		assert n_d1 + n_d2 == len(moleculeSet)

		# Randomly boolean index molecules in mother such that each daugher gets amount calculated above
		d1_bool = np.zeros(len(moleculeSet), dtype = bool)
		d2_bool = np.zeros(len(moleculeSet), dtype = bool)
		if chromosomeToDaughter1:
			d1_bool[:] = True
		else:
			d1_bool[:] = False
		d2_bool = np.logical_not(d1_bool)

		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
			d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

		d1_unique_molecules_container.objectsNew('dnaPolymerase', n_d1, **d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('dnaPolymerase', n_d2, **d2_dividedAttributesDict)

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
