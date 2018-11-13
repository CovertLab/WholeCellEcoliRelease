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

from wholecell.utils.constants import SERIALIZED_INHERITED_STATE
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
	Divide simulated states (chromosome, bulkMolecules, and uniqueMolecules)
	of a dividing cell randomly into two daughter cells, saving the data for
	the daughter cells' `_initialConditionsFunction()` method to read.
	"""
	# Assign data from simulation required
	randomState = sim.randomState

	bulkMolecules = sim.internal_states['BulkMolecules']
	uniqueMolecules = sim.internal_states['UniqueMolecules']

	sim_data = sim.get_sim_data()

	# TODO (Eran): division should be based on both nutrient and gene perturbation condition
	current_nutrients = sim.external_states['Environment'].nutrients

	# Create output directories
	d1_path = filepath.makedirs(sim._outputDir, "Daughter1")
	d2_path = filepath.makedirs(sim._outputDir, "Daughter2")
	print('Writing daughter cell data to {} et al.'.format(d1_path))
	# TODO(jerry): In a multi-scale sim, set inherited_state_path=d1_path and
	# inherited_state_path=d2_path in the daughter cell agent_config dicts,
	# along with the correct variant_type, variant_index, and seed.

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
			bulkMolecules, uniqueMolecules, randomState, chromosome_counts,
			sim_data)
		d1_uniqueMolCntr, d2_uniqueMolCntr, daughter_elng_rates = (
			divideUniqueMolecules(uniqueMolecules, randomState,
			chromosome_counts, current_nutrients, sim)
			)

	# Save the daughter initialization state.
	initial_time = sim.time() + sim.timeStepSec()
	save_inherited_state(
		d1_path,
		is_dead=isDead,
		initial_time=initial_time,
		elng_rate=daughter_elng_rates["d1_elng_rate"],
		elng_rate_factor=daughter_elng_rates["d1_elng_rate_factor"],
		bulk_molecules=d1_bulkMolCntr,
		unique_molecules=d1_uniqueMolCntr,
		)
	save_inherited_state(
		d2_path,
		is_dead=isDead,
		initial_time=initial_time,
		elng_rate=daughter_elng_rates["d2_elng_rate"],
		elng_rate_factor=daughter_elng_rates["d2_elng_rate_factor"],
		bulk_molecules=d2_bulkMolCntr,
		unique_molecules=d2_uniqueMolCntr,
		)

	return [d1_path, d2_path]


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

	d1_chromosome_mask = np.zeros(full_chromosome_count, dtype=np.bool)
	d1_chromosome_mask[d1_chromosome_indexes] = 1

	return {"d1_chromosome_count": d1_chromosome_count,
		"d2_chromosome_count": d2_chromosome_count,
		"d1_chromosome_indexes": d1_chromosome_indexes,
		"d1_chromosome_mask": d1_chromosome_mask,
		}


def divideBulkMolecules(bulkMolecules, uniqueMolecules, randomState,
		chromosome_counts, sim_data):
	"""
	Divide bulk molecules into two daughter cells based on the division ID of
	the molecule.
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

	# Divide full chromosomes
	d1_chromosome_count = chromosome_counts['d1_chromosome_count']
	d2_chromosome_count = chromosome_counts['d2_chromosome_count']
	d1_bulk_molecules_container.countIs(d1_chromosome_count,
		bulkMolecules.divisionIds['fullChromosome'][0])
	d2_bulk_molecules_container.countIs(d2_chromosome_count,
		bulkMolecules.divisionIds['fullChromosome'][0])

	# Calculate gene copy numbers for each chromosome from chromosome state
	replicationCoordinate = sim_data.process.transcription.rnaData[
		"replicationCoordinate"]
	replication_forks = uniqueMolecules.container.objectsInCollection(
		'dnaPolymerase')
	d1_chromosome_mask = chromosome_counts['d1_chromosome_mask']

	full_chromosome_count = len(d1_chromosome_mask)

	# Initialize all gene counts to one
	gene_counts_array = np.ones(
		(full_chromosome_count, len(replicationCoordinate)),
		dtype=np.int
		)

	# Add to gene counts if there are active replication forks
	if len(replication_forks) > 0:
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = replication_forks.attrs(
			'sequenceIdx', 'sequenceLength', 'replicationRound', 'chromosomeIndex'
			)

		for chrom_idx in np.arange(0, full_chromosome_count):
			chromosomeMatch = (chromosomeIndex == chrom_idx)

			forward_fork_coordinates = sequenceLength[
				np.logical_and(chromosomeMatch, (sequenceIdx == 0))
				]
			reverse_fork_coordinates = np.negative(sequenceLength[
				np.logical_and(chromosomeMatch, (sequenceIdx == 1))
				])

			assert len(forward_fork_coordinates) == len(reverse_fork_coordinates)

			for (forward, reverse) in izip(forward_fork_coordinates,
					reverse_fork_coordinates):
				gene_counts_array[chrom_idx,
					np.logical_and(replicationCoordinate < forward,
						replicationCoordinate > reverse)
					] += 1

	# Divide genes counts based on chromosome division results
	molecule_counts = bulkMolecules.container.counts(
		bulkMolecules.divisionIds['geneCopyNumber'])

	assert np.all(molecule_counts == gene_counts_array.sum(axis=0))

	d1_gene_counts = gene_counts_array[d1_chromosome_mask, :].sum(axis=0)
	d2_gene_counts = gene_counts_array[~d1_chromosome_mask, :].sum(axis=0)

	# Set gene copy numbers in daughter cells
	d1_bulk_molecules_container.countsIs(d1_gene_counts,
		bulkMolecules.divisionIds['geneCopyNumber'])
	d2_bulk_molecules_container.countsIs(d2_gene_counts,
		bulkMolecules.divisionIds['geneCopyNumber'])

	# Divide bound TFs based on gene copy numbers for each chromosome
	molecule_counts = bulkMolecules.container.counts(
		bulkMolecules.divisionIds['boundTF'])

	# Bound TFs are first divided binomially. Since the number of bound TFs
	# cannot be larger than the copy number of the target gene, the "overshoot"
	# counts for each daughter are added to the counts for the other daughter.
	d1_tf_counts = randomState.binomial(molecule_counts, p=BINOMIAL_COEFF)
	d2_tf_counts = molecule_counts - d1_tf_counts

	rna_index_to_bound_tf_mapping = sim_data.process.transcription_regulation.rna_index_to_bound_tf_mapping

	d1_tf_overshoot = (
		d1_tf_counts - d1_gene_counts[rna_index_to_bound_tf_mapping]
		).clip(min=0)
	d2_tf_overshoot = (
		d2_tf_counts - d2_gene_counts[rna_index_to_bound_tf_mapping]
		).clip(min=0)

	d1_tf_counts = d1_tf_counts - d1_tf_overshoot + d2_tf_overshoot
	d2_tf_counts = d2_tf_counts - d2_tf_overshoot + d1_tf_overshoot

	assert np.all(d1_tf_counts + d2_tf_counts == molecule_counts)

	# Set bound TF counts in daughter cells
	d1_bulk_molecules_container.countsIs(d1_tf_counts,
		bulkMolecules.divisionIds['boundTF'])
	d2_bulk_molecules_container.countsIs(d2_tf_counts,
		bulkMolecules.divisionIds['boundTF'])

	return d1_bulk_molecules_container, d2_bulk_molecules_container


def divideUniqueMolecules(uniqueMolecules, randomState, chromosome_counts,
		current_nutrients, sim):
	"""
	Divides unique molecules of the mother cell to the two daughter cells. Each
	class of unique molecules is divided in a different way.
	- active RNA polymerases: random binomial division
	- active ribosome: random binomial division, but the ribosome elongation
	rates of the daughter cells are set such that the two daughter cells have
	equal translational capacities, with an optional noise.
	- active DNA polymerases, replisomes, oriCs: divided based on the
	chromosome each molecule is associated to.
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
		'activeReplisome', 'fullChromosome', 'activeRibosome']

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

	# Divide active replisomes according to chromosomes they were bound to
	moleculeSet = uniqueMolecules.container.objectsInCollection('activeReplisome')
	moleculeAttributeDict = uniqueMoleculesToDivide['activeReplisome']
	n_replisomes = len(moleculeSet)

	if n_replisomes > 0:
		replicationRound, chromosomeIndex = moleculeSet.attrs(
			'replicationRound', 'chromosomeIndex'
		)

		# Divide replisomes based on their chromosome indexes
		d1_bool = np.zeros(n_replisomes, dtype=bool)
		for index in d1_chromosome_indexes:
			d1_bool = np.logical_or(d1_bool, chromosomeIndex == index)
		d2_bool = np.logical_not(d1_bool)

		# Add the divided replisomes to the daughter cell containers
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

		d1_unique_molecules_container.objectsNew('activeReplisome', n_d1,
			**d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('activeReplisome', n_d2,
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


def save_inherited_state(daughter_path, **inherited_state):
	"""Save the `inherited_state` dict for a daughter cell."""
	with open(os.path.join(daughter_path, SERIALIZED_INHERITED_STATE), 'wb') as f:
		cPickle.dump(inherited_state, f, cPickle.HIGHEST_PROTOCOL)

def load_inherited_state(daughter_path):
	"""Load `inherited_state` dict to initialize a daughter cell."""
	with open(os.path.join(daughter_path, SERIALIZED_INHERITED_STATE), "rb") as f:
		inherited_state = cPickle.load(f)
	return inherited_state


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
