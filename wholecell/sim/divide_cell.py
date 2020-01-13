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
	current_media_id = sim.external_states['Environment'].current_media_id

	# Create the output directory
	sim_out_dir = filepath.makedirs(sim._outputDir)
	d1_path = os.path.join(sim_out_dir, SERIALIZED_INHERITED_STATE % 1)
	d2_path = os.path.join(sim_out_dir, SERIALIZED_INHERITED_STATE % 2)
	print('Writing daughter cell data to {} et al.'.format(d1_path))

	# Check if the cell is dead
	isDead = False
	if uniqueMolecules.container.counts(["fullChromosome"])[0] == 0 and (
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
		d1_uniqueMolCntr, d2_uniqueMolCntr, daughter_elng_rates = (
			divideUniqueMolecules(uniqueMolecules, randomState,
				chromosome_division_results, current_media_id, sim))

	# Save the daughter initialization state.
	# TODO(jerry): Include the variant_type and variant_index? The seed?
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
	full_chromosomes = uniqueMolecules.container.objectsInCollection("fullChromosome")
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
	assert d1_all_domain_indexes.size + d2_all_domain_indexes.size == domain_index_domains.size

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


def divideUniqueMolecules(uniqueMolecules, randomState, chromosome_division_results,
		current_media_id, sim):
	"""
	Divides unique molecules of the mother cell to the two daughter cells.
	There are currently two different "division modes" by which unique
	molecules can be divided.

	- active_ribosome: random binomial division, but the ribosome elongation
	rates of the daughter cells are set such that the two daughter cells have
	equal translational capacities, with an optional noise.
	- domain_index: molecules are divided based on the chromosome domain each
	molecule is associated to.
	"""

	# Initialize containers for daughter cells
	d1_unique_molecules_container = uniqueMolecules.container.emptyLike()
	d2_unique_molecules_container = uniqueMolecules.container.emptyLike()

	uniqueMoleculesToDivide = deepcopy(uniqueMolecules.uniqueMoleculeDefinitions)

	# Get indexes of chromosome domains assigned to each daughter
	d1_all_domain_indexes = chromosome_division_results['d1_all_domain_indexes']
	d2_all_domain_indexes = chromosome_division_results['d2_all_domain_indexes']

	for molecule_name, molecule_attribute_dict in uniqueMoleculesToDivide.iteritems():

		molecule_set = uniqueMolecules.container.objectsInCollection(
			molecule_name)
		molecule_attribute_dict = uniqueMoleculesToDivide[molecule_name]
		n_molecules = len(molecule_set)

		if molecule_name in uniqueMolecules.division_mode['active_ribosome']:
			# Binomially divide active ribosomes, but also set the ribosome
			# elongation rates of daughter cells such that the two daughters
			# have identical translational capacities.
			daughter_elng_rates = zero_elongation_rate()

			if n_molecules > 0:

				# Read the expected ribosome elongation rate for this environment
				sim_data = sim.get_sim_data()
				elngRate = np.min([sim_data.process.translation.ribosomeElongationRateDict[
					current_media_id].asNumber(units.aa / units.s), 21.])

				# If growth rate noise is set to True, multiply noise parameter
				# to translation capacity
				noiseMultiplier = 1.
				if sim._growthRateNoise:
					noiseMultiplier = randomState.normal(1, 0.25)

				# Calculate total translation capacity of the mother cell
				translationCapacity = elngRate * n_molecules * noiseMultiplier

				# Binomially split the total count of active ribosomes
				n_d1 = randomState.binomial(n_molecules, p=BINOMIAL_COEFF)
				n_d2 = n_molecules - n_d1
				assert n_d1 + n_d2 == n_molecules

				# Calculate ribosome elongation rates of daughter cells,
				# assuming the translation capacity was split evenly
				d1_rib_elng_rate = np.min([(translationCapacity / 2) / n_d1, 21.])
				d2_rib_elng_rate = np.min([(translationCapacity / 2) / n_d2, 21.])

				daughter_elng_rates = {
					"d1_elng_rate": d1_rib_elng_rate,
					"d2_elng_rate": d2_rib_elng_rate,
					"d1_elng_rate_factor": noiseMultiplier,
					"d2_elng_rate_factor": noiseMultiplier,
					}

				# Randomly index molecules in the mother cell with a boolean
				# value such that each daughter gets amount calculated above
				d1_bool = np.zeros(n_molecules, dtype=bool)
				d1_indexes = randomState.choice(
					range(n_molecules), size=n_d1, replace=False
					)
				d1_bool[d1_indexes] = True
				d2_bool = np.logical_not(d1_bool)
			else:
				continue

		elif molecule_name in uniqueMolecules.division_mode['domain_index']:
			# Divide molecules associated with chromosomes based on the index
			# of the chromosome domains the molecules are associated with
			if n_molecules > 0:
				domain_index = molecule_set.attr("domain_index")

				# Divide molecule based on their domain indexes
				d1_bool = np.isin(domain_index, d1_all_domain_indexes)
				d2_bool = np.isin(domain_index, d2_all_domain_indexes)

				n_d1 = d1_bool.sum()
				n_d2 = d2_bool.sum()

				assert n_molecules == n_d1 + n_d2
			else:
				continue

		else:
			raise Exception, "Division mode not specified for unique molecule %s. Unable to divide cell." % (molecule_name, )

		assert n_d1 + n_d2 == n_molecules

		# Add the divided unique molecules to the daughter cell containers
		d1_divided_attributes_dict = {}
		d2_divided_attributes_dict = {}

		for molecule_attribute in molecule_attribute_dict.iterkeys():
			d1_divided_attributes_dict[molecule_attribute] = (
				molecule_set.attr(molecule_attribute)[d1_bool]
			)
			d2_divided_attributes_dict[molecule_attribute] = (
				molecule_set.attr(molecule_attribute)[d2_bool]
			)

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
