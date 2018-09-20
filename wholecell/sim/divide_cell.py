"""
divide_cell.py
Functions needed for division from mother cell into daughter cells
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import os
import cPickle
from copy import deepcopy
from itertools import izip

from wholecell.io.tablewriter import TableWriter
from wholecell.utils import filepath
from wholecell.utils import units

BINOMIAL_COEFF = 0.5


def zero_elongation_rate():
	return {
		"d1_elng_rate": 0.,
		"d2_elng_rate": 0.,
		"d1_elng_rate_factor": 0.,
		"d2_elng_rate_factor": 0.,
		}


def divide_cell(sim):
	"""
	Divides simulated states (chromosome, bulkMolecules, and uniqueMolecules)
	of a dividing cell randomly into two daughter cells.
	"""
	# Assign data from simulation required
	randomState = sim.randomState

	bulkMolecules = sim.internal_states['BulkMolecules']
	uniqueMolecules = sim.internal_states['UniqueMolecules']

	# TODO (Eran): division should be based on both nutrient and gene perturbation condition
	current_nutrients = sim.external_states['Environment'].nutrients

	# Create output directories
	d1_path = filepath.makedirs(sim._outputDir, "Daughter1")
	d2_path = filepath.makedirs(sim._outputDir, "Daughter2")
	print('Writing daughter cell data to {} et al.'.format(d1_path))
	# TODO(jerry): In a multi-scale sim, set inherited_state_path=d1_path and
	# inherited_state_path=d2_path in the daughter cell agent_config dicts,
	# along with the correct variant_type, variant_index, and seed.

	# Check for uneven numbers of partial chromosomes. This should not happen
	# too often if the four partial chromosomes are elongated in a roughly
	# synchronized way.
	# TODO (Gwanggyu): try to handle this case instead of raising an exception
	partial_chromosome_counts = bulkMolecules.container.counts(
		bulkMolecules.divisionIds['partialChromosome'])
	uneven_counts = partial_chromosome_counts - partial_chromosome_counts.min()
	if uneven_counts.any():
		raise Exception("You won the lottery! There is an uneven number of partial chromosomes...")

	# Transform any leftover partial chromosomes into full a chromosome. This
	# should have happened in the chromosome_formation process but we could get
	# unlucky and miss this in the final timestep.
	bulkMolecules.container.countInc(
		partial_chromosome_counts.min(),
		bulkMolecules.divisionIds['fullChromosome'][0]
		)

	# Check if the cell is dead
	isDead = False
	if bulkMolecules.container.count(
			bulkMolecules.divisionIds['fullChromosome'][0]) == 0 and (
			sim.time() - sim.initialTime()) > sim.lengthSec():
		# If the cell does not have any full chromosomes at the end of its
		# maximal simulation duration, the cell is considered dead
		isDead = True
	elif sim._isDead:
		isDead = True

	with open(os.path.join(sim._outputDir, "Daughter1", "IsDead.cPickle"), 'wb') as f:
		cPickle.dump(isDead, f)
	with open(os.path.join(sim._outputDir, "Daughter2", "IsDead.cPickle"), 'wb') as f:
		cPickle.dump(isDead, f)

	if isDead:
		# Cell is dead - set daughter cell containers to empty values
		d1_bulkMolCntr = bulkMolecules.container.emptyLike()
		d2_bulkMolCntr = bulkMolecules.container.emptyLike()
		d1_uniqueMolCntr = uniqueMolecules.container.emptyLike()
		d2_uniqueMolCntr = uniqueMolecules.container.emptyLike()
		daughter_elng_rates = zero_elongation_rate()
	else:
		# Divide the chromosome into two daughter cells
		# The output is used when dividing both bulk molecules and unique
		# molecules
		chromosome_counts = chromosomeDivision(bulkMolecules, randomState)

		# Create divided containers
		d1_bulkMolCntr, d2_bulkMolCntr = divideBulkMolecules(
			bulkMolecules, randomState, chromosome_counts)
		d1_uniqueMolCntr, d2_uniqueMolCntr, daughter_elng_rates = (
			divideUniqueMolecules(uniqueMolecules, randomState,
			chromosome_counts, current_nutrients, sim)
			)

	# Save divided containers
	saveContainer(d1_bulkMolCntr, os.path.join(
		sim._outputDir, "Daughter1", "BulkMolecules"))
	saveContainer(d2_bulkMolCntr, os.path.join(
		sim._outputDir, "Daughter2", "BulkMolecules"))
	saveContainer(d1_uniqueMolCntr, os.path.join(
		sim._outputDir, "Daughter1", "UniqueMolecules"))
	saveContainer(d2_uniqueMolCntr, os.path.join(
		sim._outputDir, "Daughter2", "UniqueMolecules"))

	with open(os.path.join(sim._outputDir, "Daughter1", "ElngRate.cPickle"), 'wb') as f:
		cPickle.dump(daughter_elng_rates["d1_elng_rate"], f)
	with open(os.path.join(sim._outputDir, "Daughter2", "ElngRate.cPickle"), 'wb') as f:
		cPickle.dump(daughter_elng_rates["d2_elng_rate"], f)
	with open(os.path.join(sim._outputDir, "Daughter1", "elng_rate_factor.cPickle"), 'wb') as f:
		cPickle.dump(daughter_elng_rates["d1_elng_rate_factor"], f)
	with open(os.path.join(sim._outputDir, "Daughter2", "elng_rate_factor.cPickle"), 'wb') as f:
		cPickle.dump(daughter_elng_rates["d2_elng_rate_factor"], f)

	# Save daughter cell initial time steps
	saveTime(sim.time(), os.path.join(sim._outputDir, "Daughter1", "Time"),
		sim.timeStepSec())
	saveTime(sim.time(), os.path.join(sim._outputDir, "Daughter2", "Time"),
		sim.timeStepSec())


def chromosomeDivision(bulkMolecules, randomState):
	"""
	Randomly divides full chromosome counts into two daughter cells. If there
	are even number of chromosomes, the count is divided evenly to the two
	daughter cells. If there are odd number of chromosomes, one daughter cell
	is chosen randomly to receive one more.
	Note: the actual allocation of chromosomes to the daughter cells are not
	performed in this function.
	"""
	full_chromosome_count = bulkMolecules.container.count(
		bulkMolecules.divisionIds['fullChromosome'][0])

	# Divide evenly between both daughters if there are even number of
	# chromosomes or give an extra to one cell if odd
	if full_chromosome_count % 2 == 0:
		d1_chromosome_count = full_chromosome_count // 2
		d2_chromosome_count = full_chromosome_count // 2
	else:
		d1_chromosome_count = full_chromosome_count // 2
		d1_chromosome_count += randomState.binomial(1, p=BINOMIAL_COEFF)
		d2_chromosome_count = full_chromosome_count - d1_chromosome_count

	assert d1_chromosome_count + d2_chromosome_count == full_chromosome_count

	d1_chromosome_indexes = randomState.choice(
		np.arange(full_chromosome_count),
		size=d1_chromosome_count, replace=False
		)

	return {"d1_chromosome_count": d1_chromosome_count,
		"d2_chromosome_count": d2_chromosome_count,
		"d1_chromosome_indexes": d1_chromosome_indexes
		}


def divideBulkMolecules(bulkMolecules, randomState, chromosome_counts):
	"""
	Randomly divides bulk molecules into two daughter cells based on the
	division ID of the molecule. 
	"""
	# Initialize containers for daughter cells
	d1_bulk_molecules_container = bulkMolecules.container.emptyLike()
	d2_bulk_molecules_container = bulkMolecules.container.emptyLike()

	# Binomially divide bulk molecules with division IDs set as 'binomial'
	molecule_counts = bulkMolecules.container.counts(
		bulkMolecules.divisionIds['binomial'])
	d1_counts = randomState.binomial(molecule_counts, p=BINOMIAL_COEFF)
	d2_counts = molecule_counts - d1_counts
	assert np.all(d1_counts + d2_counts == molecule_counts)

	# Set bulk molecule counts in daughter cells
	d1_bulk_molecules_container.countsIs(d1_counts,
		bulkMolecules.divisionIds['binomial'])
	d2_bulk_molecules_container.countsIs(d2_counts,
		bulkMolecules.divisionIds['binomial'])

	# Check that no bulk molecules need to be divided equally
	assert len(bulkMolecules.divisionIds['equally']) == 0

	# Handle chromosome division
	d1_chromosome_count = chromosome_counts['d1_chromosome_count']
	d2_chromosome_count = chromosome_counts['d2_chromosome_count']
	d1_bulk_molecules_container.countIs(d1_chromosome_count,
		bulkMolecules.divisionIds['fullChromosome'][0])
	d2_bulk_molecules_container.countIs(d2_chromosome_count,
		bulkMolecules.divisionIds['fullChromosome'][0])

	# Set basal transcription values to 1 for both daughter cells
	# TODO: Add assertion that these elements have zero mass
	d1_bulk_molecules_container.countsIs(1,
		bulkMolecules.divisionIds['setTo1'])
	d2_bulk_molecules_container.countsIs(1,
		bulkMolecules.divisionIds['setTo1'])

	return d1_bulk_molecules_container, d2_bulk_molecules_container


def divideUniqueMolecules(uniqueMolecules, randomState, chromosome_counts, current_nutrients, sim):
	"""
	Divides unique molecules of the mother cell to the two daughter cells. Each
	class of unique molecules is divided in a different way.
	- active RNA polymerases: random binomial division
	- active ribosome: random binomial division, but the ribosome elongation
	rates of the daughter cells are set such that the two daughter cells have
	equal translational capacities, with an optional noise.
	- active DNA polymerases, oriCs: divided based on the chromosome each
	molecule is associated to.
	- full chromosomes: assign to both daughter cells the same count of unique
	fullChromosome states as the mother cell. (Note: this unique state is just
	a placeholder for the cell division time attribute - its count essentially
	does not represent anything.)
	"""
	# Initialize containers for daughter cells
	d1_unique_molecules_container = uniqueMolecules.container.emptyLike()
	d2_unique_molecules_container = uniqueMolecules.container.emptyLike()

	uniqueMoleculesToDivide = deepcopy(uniqueMolecules.uniqueMoleculeDefinitions)

	# Get counts and indexes of chromosomes that were assigned to each daughter
	d1_chromosome_count = chromosome_counts['d1_chromosome_count']
	d2_chromosome_count = chromosome_counts['d2_chromosome_count']
	d1_chromosome_indexes = chromosome_counts['d1_chromosome_indexes']

	# List of unique molecules that should not be binomially divided
	# Note: the only unique molecule currently not in this list is active RNA
	# polymerase.
	# TODO (Gwanggyu): RNAPs should also be unequally divided once gene dosage
	# is modeled and RNAPs are assigned to a specific chromosome.
	nonbinomial_unique_molecules = ['dnaPolymerase', 'originOfReplication',
		'fullChromosome', 'activeRibosome']

	# Binomially divide unique molecules that should be binomially split
	# Note: again, the only unique molecules split here are the active RNAPs.
	for moleculeName, moleculeAttributeDict in uniqueMoleculesToDivide.iteritems():
		if moleculeName in nonbinomial_unique_molecules:
			continue

		# Get set of molecules to divide and calculate counts going to each
		# daughter cell
		moleculeSet = uniqueMolecules.container.objectsInCollection(moleculeName)

		if len(moleculeSet) > 0:
			n_d1 = randomState.binomial(len(moleculeSet), p=BINOMIAL_COEFF)
			n_d2 = len(moleculeSet) - n_d1
			assert n_d1 + n_d2 == len(moleculeSet)

			# Randomly index molecules in the mother cell with a boolean value
			# such that each daughter gets amount calculated above
			d1_bool = np.zeros(len(moleculeSet), dtype=bool)
			d1_indexes = randomState.choice(range(len(moleculeSet)), size=n_d1,
				replace=False)
			d1_bool[d1_indexes] = True
			d2_bool = np.logical_not(d1_bool)

			# Add the divided unique molecules to the daughter cell containers
			d1_dividedAttributesDict = {}
			d2_dividedAttributesDict = {}
			for moleculeAttribute in moleculeAttributeDict.iterkeys():
				d1_dividedAttributesDict[moleculeAttribute] = (
					moleculeSet.attr(moleculeAttribute)[d1_bool]
				)
				d2_dividedAttributesDict[moleculeAttribute] = (
					moleculeSet.attr(moleculeAttribute)[d2_bool]
				)

			d1_unique_molecules_container.objectsNew(moleculeName, n_d1,
				**d1_dividedAttributesDict)
			d2_unique_molecules_container.objectsNew(moleculeName, n_d2,
				**d2_dividedAttributesDict)

	# Binomially divide active ribosomes, but also set the ribosome elongation
	# rates of daughter cells such that the two daughters have identical
	# translation capacities.
	moleculeSet = uniqueMolecules.container.objectsInCollection('activeRibosome')
	moleculeAttributeDict = uniqueMoleculesToDivide['activeRibosome']
	n_ribosomes = len(moleculeSet)

	daughter_elng_rates = zero_elongation_rate()

	if n_ribosomes > 0:
		# Read the ribosome elongation rate of the mother cell
		polypeptide_elongation = sim.processes["PolypeptideElongation"]
		elngRate = np.min([polypeptide_elongation.ribosomeElongationRateDict[
			current_nutrients].asNumber(units.aa / units.s), 21.])

		# If growth rate noise is set to True, multiply noise parameter to
		# translation capacity
		noiseMultiplier = 1.
		if sim._growthRateNoise:
			noiseMultiplier = randomState.normal(1, 0.25)

		# Calculate total translation capacity of the mother cell
		translationCapacity = elngRate * n_ribosomes * noiseMultiplier

		# Binomially split the total count of active ribosomes
		n_d1 = randomState.binomial(n_ribosomes, p=BINOMIAL_COEFF)
		n_d2 = n_ribosomes - n_d1
		assert n_d1 + n_d2 == n_ribosomes

		# Calculate ribosome elongation rates of daughter cells, assuming the
		# translation capacity was split evenly
		d1_rib_elng_rate = np.min([(translationCapacity / 2) / n_d1, 21.])
		d2_rib_elng_rate = np.min([(translationCapacity / 2) / n_d2, 21.])

		daughter_elng_rates = {
			"d1_elng_rate": d1_rib_elng_rate,
			"d2_elng_rate": d2_rib_elng_rate,
			"d1_elng_rate_factor": noiseMultiplier,
			"d2_elng_rate_factor": noiseMultiplier,
			}

		# Randomly index molecules in the mother cell with a boolean value
		# such that each daughter gets amount calculated above
		d1_bool = np.zeros(n_ribosomes, dtype=bool)
		d1_indexes = randomState.choice(range(n_ribosomes), size=n_d1, replace=False)
		d1_bool[d1_indexes] = True
		d2_bool = np.logical_not(d1_bool)

		# Add the divided ribosomes to the daughter cell containers
		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = (
				moleculeSet.attr(moleculeAttribute)[d1_bool]
			)
			d2_dividedAttributesDict[moleculeAttribute] = (
				moleculeSet.attr(moleculeAttribute)[d2_bool]
			)

		d1_unique_molecules_container.objectsNew('activeRibosome', n_d1,
			**d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('activeRibosome', n_d2,
			**d2_dividedAttributesDict)

	# Divide DNA polymerases according to chromosomes they were bound to
	moleculeSet = uniqueMolecules.container.objectsInCollection('dnaPolymerase')
	moleculeAttributeDict = uniqueMoleculesToDivide['dnaPolymerase']
	n_dnaps = len(moleculeSet)

	if n_dnaps > 0:
		sequenceIdx, replicationRound, chromosomeIndex = moleculeSet.attrs(
			'sequenceIdx', 'replicationRound', 'chromosomeIndex'
		)

		# Divide polymerases based on their chromosome indexes
		d1_bool = np.zeros(n_dnaps, dtype=bool)
		for index in d1_chromosome_indexes:
			d1_bool = np.logical_or(d1_bool, chromosomeIndex == index)
		d2_bool = np.logical_not(d1_bool)

		# Add the divided DNAPs to the daughter cell containers
		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = (
				moleculeSet.attr(moleculeAttribute)[d1_bool]
			)
			d2_dividedAttributesDict[moleculeAttribute] = (
				moleculeSet.attr(moleculeAttribute)[d2_bool]
			)

		n_d1 = d1_bool.sum()
		n_d2 = d2_bool.sum()

		# Reset the chromosome indexes of the polymerases assigned to each daughter cell
		d1_dividedAttributesDict['chromosomeIndex'] = resetChromosomeIndex(
			d1_dividedAttributesDict['chromosomeIndex'], d1_chromosome_count)
		d2_dividedAttributesDict['chromosomeIndex'] = resetChromosomeIndex(
			d2_dividedAttributesDict['chromosomeIndex'], d2_chromosome_count)

		d1_unique_molecules_container.objectsNew('dnaPolymerase', n_d1,
			**d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('dnaPolymerase', n_d2,
			**d2_dividedAttributesDict)

	# Divide oriCs according to the chromosomes they are associated to
	moleculeSet = uniqueMolecules.container.objectsInCollection('originOfReplication')
	moleculeAttributeDict = uniqueMoleculesToDivide['originOfReplication']
	n_oric = len(moleculeSet)

	if n_oric > 0:
		chromosomeIndex = moleculeSet.attr('chromosomeIndex')

		# Divide oriC's based on their chromosome indexes
		d1_bool = np.zeros(n_oric, dtype=bool)
		for index in d1_chromosome_indexes:
			d1_bool = np.logical_or(d1_bool, chromosomeIndex == index)
		d2_bool = np.logical_not(d1_bool)

		# Add the divided OriCs to the daughter cell containers
		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = (
				moleculeSet.attr(moleculeAttribute)[d1_bool]
			)
			d2_dividedAttributesDict[moleculeAttribute] = (
				moleculeSet.attr(moleculeAttribute)[d2_bool]
			)

		n_d1 = d1_bool.sum()
		n_d2 = d2_bool.sum()

		# Reset chromosome indexes of polymerases assigned to each daughter
		d1_dividedAttributesDict['chromosomeIndex'] = resetChromosomeIndex(
			d1_dividedAttributesDict['chromosomeIndex'], d1_chromosome_count)
		d2_dividedAttributesDict['chromosomeIndex'] = resetChromosomeIndex(
			d2_dividedAttributesDict['chromosomeIndex'], d2_chromosome_count)

		d1_unique_molecules_container.objectsNew('originOfReplication', n_d1,
			**d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('originOfReplication', n_d2,
			**d2_dividedAttributesDict)

	# Divide unique fullChromosomes
	moleculeSet = uniqueMolecules.container.objectsInCollection('fullChromosome')
	moleculeAttributeDict = uniqueMoleculesToDivide['fullChromosome']
	n_full_chromosome = len(moleculeSet)

	if n_full_chromosome > 0:
		# Add full chromosome unique molecule states to both daughter cells
		# Note: The fullChromosome unique molecule is simply a placeholder to
		# log cell division data. Daughter cells should get one molecule each.
		d1_bool = np.ones(n_full_chromosome, dtype=bool)
		d2_bool = np.ones(n_full_chromosome, dtype=bool)

		n_d1 = n_full_chromosome
		n_d2 = n_full_chromosome

		d1_dividedAttributesDict = {}
		d2_dividedAttributesDict = {}
		for moleculeAttribute in moleculeAttributeDict.iterkeys():
			d1_dividedAttributesDict[moleculeAttribute] = (
				moleculeSet.attr(moleculeAttribute)[d1_bool]
			)
			d2_dividedAttributesDict[moleculeAttribute] = (
				moleculeSet.attr(moleculeAttribute)[d2_bool]
			)

		d1_unique_molecules_container.objectsNew('fullChromosome', n_d1,
			**d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('fullChromosome', n_d2,
			**d2_dividedAttributesDict)

	return d1_unique_molecules_container, d2_unique_molecules_container, daughter_elng_rates


def saveContainer(container, path):
	table_writer = TableWriter(path)
	container.tableCreate(table_writer)
	container.tableAppend(table_writer)


def saveTime(finalTime, path, timeStepSec):
	timeFile = TableWriter(path)

	# Metadata
	timeFile.writeAttributes(initialTime=finalTime + timeStepSec)
	timeFile.close()


def resetChromosomeIndex(oldChromosomeIndex, chromosomeCount):
	"""
	Resets the chromosome index array of the daughter cell such that the array
	elements index the chromosomes with consecutive integers starting from
	zero. Returns the new index array.
	"""
	newChromosomeIndex = np.zeros_like(oldChromosomeIndex, dtype=np.int)
	for newIndex, oldIndex in izip(np.arange(chromosomeCount), np.unique(oldChromosomeIndex)):
		indexMatch = (oldChromosomeIndex == oldIndex)
		newChromosomeIndex[indexMatch] = newIndex

	return newChromosomeIndex
