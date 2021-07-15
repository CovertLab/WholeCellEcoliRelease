"""
divide_cell.py
Functions needed for division from mother cell into daughter cells
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import os
from six.moves import cPickle
from copy import deepcopy

from wholecell.utils.constants import SERIALIZED_INHERITED_STATE
from wholecell.utils import filepath
from six.moves import zip

BINOMIAL_COEFF = 0.5


class UniqueMoleculeDivisionError(Exception):
	pass


def zero_elongation_rate():
	return {
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

	# Create the output directory
	sim_out_dir = filepath.makedirs(sim._outputDir)
	d1_path = os.path.join(sim_out_dir, SERIALIZED_INHERITED_STATE % 1)
	d2_path = os.path.join(sim_out_dir, SERIALIZED_INHERITED_STATE % 2)
	print('Writing daughter cell data to {} et al.'.format(d1_path))

	# Check if the cell is dead
	isDead = False
	if uniqueMolecules.container.counts(['full_chromosome'])[0] == 0 and (
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
		# Divide full chromosomes into two daughter cells
		# The output is used when dividing both bulk molecules and unique
		# molecules
		no_child_place_holder = sim_data.process.replication.no_child_place_holder
		chromosome_division_results = chromosomeDivision(
			uniqueMolecules, randomState, no_child_place_holder)

		# Create divided containers
		d1_bulkMolCntr, d2_bulkMolCntr = divideBulkMolecules(
			bulkMolecules, randomState)
		d1_uniqueMolCntr, d2_uniqueMolCntr, daughter_elng_rates = divideUniqueMolecules(
			uniqueMolecules, randomState, chromosome_division_results, sim)

	# Save the daughter initialization state.
	# TODO(jerry): Include the variant_type and variant_index? The seed?
	initial_time = sim.time() + sim.timeStepSec()
	save_inherited_state(
		d1_path,
		is_dead=isDead,
		initial_time=initial_time,
		elng_rate_factor=daughter_elng_rates["d1_elng_rate_factor"],
		bulk_molecules=d1_bulkMolCntr,
		unique_molecules=d1_uniqueMolCntr,
		)
	save_inherited_state(
		d2_path,
		is_dead=isDead,
		initial_time=initial_time,
		elng_rate_factor=daughter_elng_rates["d2_elng_rate_factor"],
		bulk_molecules=d2_bulkMolCntr,
		unique_molecules=d2_uniqueMolCntr,
		)

	return [d1_path, d2_path]


def chromosomeDivision(uniqueMolecules, randomState, no_child_place_holder):
	"""
	Splits chromosome domain indexes into two daughter cells. If there are an
	even number of full chromosomes, each cell gets an equal amount of full
	chromosomes, and all the descendent domains of the oldest domains of these
	full chromosomes. If there is an odd number of full chromosomes, one cell
	gets one more full chromosome than the other.
	Note: the actual allocation of full chromosome molecules to the daughter
	cells are not performed in this function.
	"""
	# Read attributes of full chromosomes and chromosome domains
	full_chromosomes = uniqueMolecules.container.objectsInCollection('full_chromosome')
	domain_index_full_chroms = full_chromosomes.attr("domain_index")

	chromosome_domains = uniqueMolecules.container.objectsInCollection("chromosome_domain")
	domain_index_domains, child_domains = chromosome_domains.attrs(
		"domain_index", "child_domains"
		)

	# Randomly decide which daughter gets first full chromosome
	d1_gets_first_chromosome = randomState.rand() < BINOMIAL_COEFF

	index = not d1_gets_first_chromosome
	d1_domain_index_full_chroms = domain_index_full_chroms[index::2]
	d2_domain_index_full_chroms = domain_index_full_chroms[not index::2]
	d1_all_domain_indexes = get_descendent_domains(
		d1_domain_index_full_chroms, domain_index_domains,
		child_domains, no_child_place_holder
		)
	d2_all_domain_indexes = get_descendent_domains(
		d2_domain_index_full_chroms, domain_index_domains,
		child_domains, no_child_place_holder
		)

	# Check that the domains are being divided correctly
	assert np.intersect1d(d1_all_domain_indexes, d2_all_domain_indexes).size == 0

	d1_chromosome_count = d1_domain_index_full_chroms.size
	d2_chromosome_count = d2_domain_index_full_chroms.size

	return {
		"d1_all_domain_indexes": d1_all_domain_indexes,
		"d2_all_domain_indexes": d2_all_domain_indexes,
		"d1_chromosome_count": d1_chromosome_count,
		"d2_chromosome_count": d2_chromosome_count,
		}


def divideBulkMolecules(bulkMolecules, randomState):
	"""
	Divide bulk molecules into two daughter cells based on the division mode of
	the molecule.
	"""
	# Initialize containers for daughter cells
	d1_bulk_molecules_container = bulkMolecules.container.emptyLike()
	d2_bulk_molecules_container = bulkMolecules.container.emptyLike()

	# Binomially divide bulk molecules with division mode set as 'binomial'
	molecule_counts = bulkMolecules.container.counts(
		bulkMolecules.division_mode['binomial'])
	d1_counts = randomState.binomial(molecule_counts, p=BINOMIAL_COEFF)
	d2_counts = molecule_counts - d1_counts
	assert np.all(d1_counts + d2_counts == molecule_counts)

	# Set bulk molecule counts in daughter cells
	d1_bulk_molecules_container.countsIs(d1_counts,
		bulkMolecules.division_mode['binomial'])
	d2_bulk_molecules_container.countsIs(d2_counts,
		bulkMolecules.division_mode['binomial'])

	# Check that no bulk molecules need to be divided equally
	assert len(bulkMolecules.division_mode['equally']) == 0

	return d1_bulk_molecules_container, d2_bulk_molecules_container


def divideUniqueMolecules(uniqueMolecules, randomState,
		chromosome_division_results, sim):
	"""
	Divides unique molecules of the mother cell to the two daughter cells.
	There are currently three different "division modes" by which unique
	molecules can be divided.

	- domain_index: Molecules are divided based on the chromosome domain each
		molecule is associated to.
	- RNA: If RNA is a full mRNA transcript, divide binomially. If RNA
		is a partial transcript, follow the chromosome domain that the
		associated RNA polymerase molecule is bound to.
	- active_ribosome: Follow the mRNA molecule that the ribosome is bound to.
		If the mRNA molecule has already been degraded, divide binomially.

	Since RNA division is dependent on RNA polymerase division, and
	active_ribosome division is dependent on RNA division, division must occur
	in the order shown above.
	"""

	# Initialize containers for daughter cells
	d1_unique_molecules_container = uniqueMolecules.container.emptyLike()
	d2_unique_molecules_container = uniqueMolecules.container.emptyLike()

	uniqueMoleculesToDivide = deepcopy(uniqueMolecules.unique_molecule_definitions)

	# Get indexes of chromosome domains assigned to each daughter
	d1_all_domain_indexes = chromosome_division_results['d1_all_domain_indexes']
	d2_all_domain_indexes = chromosome_division_results['d2_all_domain_indexes']

	# Initialize daughter cell ribosome elongation rates
	daughter_elng_rates = zero_elongation_rate()

	# Initialize arrays for unique indexes of chromosome-bound molecules. The
	# number -1 is used to indicate terC and is always retained.
	d1_chromosome_bound_unique_indexes = np.array([-1], dtype=np.int64)
	d2_chromosome_bound_unique_indexes = np.array([-1], dtype=np.int64)
	d1_chromosome_bound_new_unique_indexes = np.array([-1], dtype=np.int64)
	d2_chromosome_bound_new_unique_indexes = np.array([-1], dtype=np.int64)

	# Divide molecules with division mode "domain_index"
	for molecule_name in uniqueMolecules.division_mode['domain_index']:
		molecule_set = uniqueMolecules.container.objectsInCollection(
			molecule_name)
		molecule_attribute_dict = uniqueMoleculesToDivide[molecule_name]
		n_molecules = len(molecule_set)

		# Divide molecules associated with chromosomes based on the index
		# of the chromosome domains the molecules are associated with
		domain_index = molecule_set.attr("domain_index")

		# Divide molecule based on their domain indexes
		d1_bool = np.isin(domain_index, d1_all_domain_indexes)
		d2_bool = np.isin(domain_index, d2_all_domain_indexes)

		n_d1 = np.count_nonzero(d1_bool)
		n_d2 = np.count_nonzero(d2_bool)

		unique_index = molecule_set.attr("unique_index")
		d1_chromosome_bound_unique_indexes = np.append(
			d1_chromosome_bound_unique_indexes,
			unique_index[d1_bool])
		d2_chromosome_bound_unique_indexes = np.append(
			d2_chromosome_bound_unique_indexes,
			unique_index[d2_bool])

		# Only the chromosome domains that physically still exist are inherited
		# to daughter cells
		if molecule_name != 'chromosome_domain':
			assert n_molecules == n_d1 + n_d2

		assert np.count_nonzero(np.logical_and(d1_bool, d2_bool)) == 0

		# Add the divided unique molecules to the daughter cell containers
		d1_divided_attributes_dict, d2_divided_attributes_dict = get_divided_attributes(
			molecule_set, molecule_attribute_dict, d1_bool, d2_bool)

		d1_unique_indexes = d1_unique_molecules_container.objectsNew(
			molecule_name, n_d1,
			**d1_divided_attributes_dict)
		d2_unique_indexes = d2_unique_molecules_container.objectsNew(
			molecule_name, n_d2,
			**d2_divided_attributes_dict)

		# Remember the new unique indexes that the molecules were given
		d1_chromosome_bound_new_unique_indexes = np.append(
			d1_chromosome_bound_new_unique_indexes,
			d1_unique_indexes)
		d2_chromosome_bound_new_unique_indexes = np.append(
			d2_chromosome_bound_new_unique_indexes,
			d2_unique_indexes)

	for molecule_name in uniqueMolecules.division_mode['chromosomal_segment']:
		molecule_set = uniqueMolecules.container.objectsInCollection(
			molecule_name)
		molecule_attribute_dict = uniqueMoleculesToDivide[molecule_name]
		n_molecules = len(molecule_set)

		# Divide chromosomal segments based on the index of the chromosome
		# domains that the segments belong to
		domain_index = molecule_set.attr("domain_index")

		# Divide molecule based on their domain indexes
		d1_bool = np.isin(domain_index, d1_all_domain_indexes)
		d2_bool = np.isin(domain_index, d2_all_domain_indexes)

		n_d1 = np.count_nonzero(d1_bool)
		n_d2 = np.count_nonzero(d2_bool)

		assert n_molecules == n_d1 + n_d2
		assert np.count_nonzero(np.logical_and(d1_bool, d2_bool)) == 0

		# Add the divided unique molecules to the daughter cell containers
		d1_divided_attributes_dict, d2_divided_attributes_dict = get_divided_attributes(
			molecule_set, molecule_attribute_dict, d1_bool, d2_bool)

		# Reset the boundary_molecule_indexes attributes of chromosomal
		# segments with the new unique indexes given to the boundary molecules
		d1_old_boundary_molecule_indexes = d1_divided_attributes_dict['boundary_molecule_indexes']
		d2_old_boundary_molecule_indexes = d2_divided_attributes_dict['boundary_molecule_indexes']

		d1_new_boundary_molecule_indexes = remap_unique_indexes(
			d1_old_boundary_molecule_indexes,
			d1_chromosome_bound_unique_indexes,
			d1_chromosome_bound_new_unique_indexes)
		d2_new_boundary_molecule_indexes = remap_unique_indexes(
			d2_old_boundary_molecule_indexes,
			d2_chromosome_bound_unique_indexes,
			d2_chromosome_bound_new_unique_indexes)

		d1_divided_attributes_dict['boundary_molecule_indexes'] = d1_new_boundary_molecule_indexes
		d2_divided_attributes_dict['boundary_molecule_indexes'] = d2_new_boundary_molecule_indexes

		d1_unique_molecules_container.objectsNew(
			molecule_name, n_d1, **d1_divided_attributes_dict)
		d2_unique_molecules_container.objectsNew(
			molecule_name, n_d2, **d2_divided_attributes_dict)

	d1_RNA_unique_indexes = None
	d2_RNA_unique_indexes = None
	d1_RNA_new_unique_indexes = None
	d2_RNA_new_unique_indexes = None

	# Divide molecules with division mode "RNA"
	for molecule_name in uniqueMolecules.division_mode['RNA']:
		molecule_set = uniqueMolecules.container.objectsInCollection(
			molecule_name)
		molecule_attribute_dict = uniqueMoleculesToDivide[molecule_name]
		n_molecules = len(molecule_set)

		# Divide full mRNA transcripts binomially, and partial transcripts
		# following the chromosome domain that the associated RNA
		# polymerase molecule is bound to.
		if n_molecules > 0:
			is_full_transcript, RNAP_index, RNA_unique_index = molecule_set.attrs(
				"is_full_transcript", "RNAP_index", "unique_index")

			d1_bool = np.zeros(n_molecules, dtype=bool)
			d2_bool = np.zeros(n_molecules, dtype=bool)

			# Divide full transcripts binomially
			full_transcript_indexes = np.where(is_full_transcript)[0]
			if len(full_transcript_indexes) > 0:
				n_full_d1 = randomState.binomial(
					np.count_nonzero(is_full_transcript), p=BINOMIAL_COEFF)
				full_d1_indexes = randomState.choice(
					full_transcript_indexes, size=n_full_d1,
					replace=False)
				full_d2_indexes = np.setdiff1d(full_transcript_indexes,
					full_d1_indexes)

				d1_bool[full_d1_indexes] = True
				d2_bool[full_d2_indexes] = True

			# Divide partial transcripts based on how their associated
			# RNAPs were divided
			partial_transcript_indexes = np.where(
				np.logical_not(is_full_transcript))[0]
			RNAP_index_partial_transcripts = RNAP_index[
				partial_transcript_indexes]

			partial_d1_indexes = partial_transcript_indexes[
				np.isin(RNAP_index_partial_transcripts, d1_chromosome_bound_unique_indexes)]
			partial_d2_indexes = partial_transcript_indexes[
				np.isin(RNAP_index_partial_transcripts, d2_chromosome_bound_unique_indexes)]

			d1_bool[partial_d1_indexes] = True
			d2_bool[partial_d2_indexes] = True

			n_d1 = np.count_nonzero(d1_bool)
			n_d2 = np.count_nonzero(d2_bool)

			if molecule_name == 'RNA':
				# Save data for future use (active ribosome division)
				d1_RNA_unique_indexes = RNA_unique_index[d1_bool]
				d2_RNA_unique_indexes = RNA_unique_index[d2_bool]
		else:
			if molecule_name == 'RNA':
				d1_RNA_unique_indexes = np.array([], dtype=np.int64)
				d2_RNA_unique_indexes = np.array([], dtype=np.int64)
			continue

		assert n_molecules == n_d1 + n_d2
		assert np.count_nonzero(np.logical_and(d1_bool, d2_bool)) == 0

		# Add the divided unique molecules to the daughter cell containers
		d1_divided_attributes_dict, d2_divided_attributes_dict = get_divided_attributes(
			molecule_set, molecule_attribute_dict, d1_bool, d2_bool)

		# Reset the RNAP_index attributes of RNAs with the new unique
		# indexes given to RNAP molecules
		d1_RNAP_index_old = d1_divided_attributes_dict['RNAP_index']
		d2_RNAP_index_old = d2_divided_attributes_dict['RNAP_index']

		d1_RNAP_index_new = remap_unique_indexes(d1_RNAP_index_old,
			d1_chromosome_bound_unique_indexes, d1_chromosome_bound_new_unique_indexes)
		d2_new_boundary_molecule_indexes = remap_unique_indexes(d2_RNAP_index_old,
			d2_chromosome_bound_unique_indexes, d2_chromosome_bound_new_unique_indexes)

		d1_divided_attributes_dict['RNAP_index'] = d1_RNAP_index_new
		d2_divided_attributes_dict['RNAP_index'] = d2_new_boundary_molecule_indexes

		d1_unique_indexes = d1_unique_molecules_container.objectsNew(
			molecule_name, n_d1,
			**d1_divided_attributes_dict)
		d2_unique_indexes = d2_unique_molecules_container.objectsNew(
			molecule_name, n_d2,
			**d2_divided_attributes_dict)

		if molecule_name == 'RNA':
			d1_RNA_new_unique_indexes = d1_unique_indexes
			d2_RNA_new_unique_indexes = d2_unique_indexes

	# Check that RNAs have been properly divided
	if any(v is None for v in [d1_RNA_unique_indexes, d2_RNA_unique_indexes,
			d1_RNA_new_unique_indexes, d2_RNA_new_unique_indexes]):
		raise UniqueMoleculeDivisionError(
			'RNAs must be divided before dividing active ribosomes.')

	for molecule_name in uniqueMolecules.division_mode['active_ribosome']:
		molecule_set = uniqueMolecules.container.objectsInCollection(
			molecule_name)
		molecule_attribute_dict = uniqueMoleculesToDivide[molecule_name]
		n_molecules = len(molecule_set)

		# Divide ribosomes following the mRNA molecule that each ribosome
		# is bound to.
		if n_molecules > 0:
			# If growth rate noise is set to True, multiply noise parameter
			# to translation capacity
			noiseMultiplier = 1.
			if sim._growthRateNoise:
				noiseMultiplier = randomState.normal(1, 0.25)

			daughter_elng_rates = {
				"d1_elng_rate_factor": noiseMultiplier,
				"d2_elng_rate_factor": noiseMultiplier,
				}

			# Divide ribosomes based on their mRNA index
			mRNA_index = molecule_set.attr("mRNA_index")
			d1_bool = np.isin(mRNA_index, d1_RNA_unique_indexes)
			d2_bool = np.isin(mRNA_index, d2_RNA_unique_indexes)

			# Binomially divide indexes of mRNAs that are degraded but still
			# has bound ribosomes. Note (ggsun): This happens because mRNA
			# degradation does not abort ongoing translation of the mRNA
			degraded_mRNA_indexes = np.unique(mRNA_index[
				np.logical_not(np.logical_or(d1_bool, d2_bool))])
			n_degraded_mRNA = len(degraded_mRNA_indexes)

			if n_degraded_mRNA > 0:
				n_degraded_mRNA_d1 = randomState.binomial(
					n_degraded_mRNA, p=BINOMIAL_COEFF)
				degraded_mRNA_indexes_d1 = randomState.choice(
					degraded_mRNA_indexes, size=n_degraded_mRNA_d1, replace=False)
				degraded_mRNA_indexes_d2 = np.setdiff1d(
					degraded_mRNA_indexes, degraded_mRNA_indexes_d1)

				# Divide "lost" ribosomes based on how these mRNAs were divided
				lost_ribosomes_d1 = np.isin(mRNA_index, degraded_mRNA_indexes_d1)
				lost_ribosomes_d2 = np.isin(mRNA_index, degraded_mRNA_indexes_d2)

				d1_bool[lost_ribosomes_d1] = True
				d2_bool[lost_ribosomes_d2] = True

			n_d1 = np.count_nonzero(d1_bool)
			n_d2 = np.count_nonzero(d2_bool)
		else:
			continue

		assert n_molecules == n_d1 + n_d2
		assert np.count_nonzero(np.logical_and(d1_bool, d2_bool)) == 0

		# Add the divided unique molecules to the daughter cell containers
		d1_divided_attributes_dict, d2_divided_attributes_dict = get_divided_attributes(
			molecule_set, molecule_attribute_dict, d1_bool, d2_bool)

		# Reset the mRNA_index attributes of ribosomes with the new unique
		# indexes given to RNA molecules
		d1_mRNA_index_old = d1_divided_attributes_dict['mRNA_index']
		d2_mRNA_index_old = d2_divided_attributes_dict['mRNA_index']

		d1_mRNA_index_new = remap_unique_indexes(d1_mRNA_index_old,
			d1_RNA_unique_indexes, d1_RNA_new_unique_indexes)
		d2_mRNA_index_new = remap_unique_indexes(d2_mRNA_index_old,
			d2_RNA_unique_indexes, d2_RNA_new_unique_indexes)

		d1_divided_attributes_dict['mRNA_index'] = d1_mRNA_index_new
		d2_divided_attributes_dict['mRNA_index'] = d2_mRNA_index_new

		d1_unique_molecules_container.objectsNew(
			molecule_name, n_d1,
			**d1_divided_attributes_dict)
		d2_unique_molecules_container.objectsNew(
			molecule_name, n_d2,
			**d2_divided_attributes_dict)

	return d1_unique_molecules_container, d2_unique_molecules_container, daughter_elng_rates


def save_inherited_state(daughter_path, **inherited_state):
	"""Save the `inherited_state` dict for a daughter cell."""
	with open(daughter_path, 'wb') as f:
		cPickle.dump(inherited_state, f, cPickle.HIGHEST_PROTOCOL)


def load_inherited_state(daughter_path):
	"""Load `inherited_state` dict to initialize a daughter cell."""
	with open(daughter_path, "rb") as f:
		inherited_state = cPickle.load(f)
	return inherited_state


def flatten(l):
	"""
	Flattens a nested list into a single list.
	"""
	return sum(l, [])


def follow_domain_tree(domain, domain_index, child_domains, place_holder):
	"""
	Recursive function that returns all the descendents of a single node in
	the domain tree, including itself.
	"""
	children_nodes = child_domains[np.where(domain_index == domain)[0][0]]

	if children_nodes[0] != place_holder:
		# If the node has children, recursively run function on each of the
		# node's two children
		branches = flatten([
			follow_domain_tree(child, domain_index, child_domains, place_holder)
			for child in children_nodes])

		# Append index of the node itself
		branches.append(domain)
		return branches

	else:
		# If the node has no children, return the index of itself
		return [domain]


def get_descendent_domains(root_domains, domain_index, child_domains, place_holder):
	"""
	Returns an array of domain indexes that are descendents of the indexes
	listed in root_domains, including the indexes in root_domains themselves.
	"""
	return np.array(flatten([
		follow_domain_tree(root_domain, domain_index, child_domains, place_holder)
		for root_domain in root_domains]))


def get_divided_attributes(
		molecule_set, molecule_attribute_dict, d1_bool, d2_bool):
	"""
	Divide the attributes of a molecule set into two separate dictionaries.
	"""
	d1_divided_attributes_dict = {}
	d2_divided_attributes_dict = {}

	for molecule_attribute in molecule_attribute_dict:
		d1_divided_attributes_dict[molecule_attribute] = (
			molecule_set.attr(molecule_attribute)[d1_bool]
		)
		d2_divided_attributes_dict[molecule_attribute] = (
			molecule_set.attr(molecule_attribute)[d2_bool]
		)

	return d1_divided_attributes_dict, d2_divided_attributes_dict


def remap_unique_indexes(indexes, old_unique_indexes, new_unique_indexes):
	"""
	Replaces the numbers in the indexes array from the numbers found in the
	old_unique_indexes array to the numbers in the new_unique_indexes array.
	If the number is not found in old_unique_indexes, the number is replaced
	with a -1. Shape of the original indexes array is retained.
	"""
	old_to_new_index = {old: new for old, new
		in zip(old_unique_indexes, new_unique_indexes)}
	new_indexes = np.array(
		[old_to_new_index.get(i, -1) for i in indexes.flatten()]
		).reshape(indexes.shape)

	return new_indexes
