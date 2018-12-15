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
		chromosome_division_results = chromosomeDivision(
			uniqueMolecules, randomState, sim_data)

		# Create divided containers
		d1_bulkMolCntr, d2_bulkMolCntr = divideBulkMolecules(
			bulkMolecules, uniqueMolecules, randomState, chromosome_division_results,
			sim_data)
		d1_uniqueMolCntr, d2_uniqueMolCntr, daughter_elng_rates = (
			divideUniqueMolecules(uniqueMolecules, randomState,
			chromosome_division_results, current_nutrients, sim)
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


def chromosomeDivision(uniqueMolecules, randomState, sim_data):
	"""
	Splits chromosome domain indexes into two daughter cells. If there are an
	even number of full chromosomes, each cell gets an equal amount of full
	chromosomes, and all the descendent domains of the mother domain of these
	full chromosomes. If there is an odd number of full chromosomes, one cell
	gets one more full chromosome than the other.
	Note: the actual allocation of full chromosome molecules to the daughter
	cells are not performed in this function.
	"""
	# Read attributes of full chromosomes and chromosome domains
	full_chromosomes = uniqueMolecules.container.objectsInCollection("fullChromosome")
	mother_domain_index = full_chromosomes.attr("mother_domain_index")
	full_chromosome_count = mother_domain_index.size

	chromosome_domains = uniqueMolecules.container.objectsInCollection("chromosome_domain")
	domain_index, child_domains = chromosome_domains.attrs(
		"domain_index", "child_domains"
		)

	# Get placeholder value for domains without children
	place_holder = sim_data.process.replication.no_child_place_holder

	# Randomly decide which daughter gets first full chromosome
	d1_gets_first_chromosome = randomState.rand() < BINOMIAL_COEFF

	if d1_gets_first_chromosome:
		d1_mother_domain_indexes = mother_domain_index[0::2]
		d1_all_domain_indexes = get_descendent_domains(
			d1_mother_domain_indexes, domain_index, child_domains, place_holder
			)
	else:
		d1_mother_domain_indexes = mother_domain_index[1::2]
		d1_all_domain_indexes = get_descendent_domains(
			d1_mother_domain_indexes, domain_index, child_domains, place_holder
			)

	d1_chromosome_count = d1_mother_domain_indexes.size
	d2_chromosome_count = full_chromosome_count - d1_chromosome_count

	return {
		"d1_mother_domain_indexes": d1_mother_domain_indexes,
		"d1_all_domain_indexes": d1_all_domain_indexes,
		"d1_chromosome_count": d1_chromosome_count,
		"d2_chromosome_count": d2_chromosome_count,
		}


def divideBulkMolecules(bulkMolecules, uniqueMolecules, randomState,
		chromosome_division_results, sim_data):
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

	# Calculate gene copy numbers for each chromosome from chromosome state
	# Get total gene copy number
	total_gene_counts = bulkMolecules.container.counts(
		bulkMolecules.divisionIds['geneCopyNumber'])

	# Get coordinates of genes and existing replication forks
	replicationCoordinate = sim_data.process.transcription.rnaData[
		"replicationCoordinate"]
	active_replisomes = uniqueMolecules.container.objectsInCollection(
		"active_replisome")

	# Get all domain indexes that should be moved to daughter 1
	d1_all_domain_indexes = chromosome_division_results['d1_all_domain_indexes']

	# Initialize gene counts for daughter 1 to one
	d1_gene_counts = np.ones(len(replicationCoordinate), dtype=np.int64)

	# Add to gene counts if there are active replication forks
	if len(active_replisomes) > 0:
		coordinates, right_replichore, domain_index = active_replisomes.attrs(
			"coordinates", "right_replichore", "domain_index"
			)

		# Get indexes of d1 domains in domain_index
		d1_domains_in_domain_index = np.where(np.in1d(domain_index, d1_all_domain_indexes))[0]

		# Filter out attributes of replisomes that belong to daughter 1
		d1_coordinates = coordinates[d1_domains_in_domain_index]
		d1_right_replichore = right_replichore[d1_domains_in_domain_index]

		forward_fork_coordinates = d1_coordinates[d1_right_replichore]
		reverse_fork_coordinates = d1_coordinates[~d1_right_replichore]

		assert len(forward_fork_coordinates) == len(reverse_fork_coordinates)

		# Increment gene counts if gene is between two replication forks
		for (forward, reverse) in izip(forward_fork_coordinates,
				reverse_fork_coordinates):
			d1_gene_counts[
				np.logical_and(replicationCoordinate < forward,
					replicationCoordinate > reverse)
				] += 1

	# Set gene copy numbers in daughter 2 by subtracting from total count
	d2_gene_counts = total_gene_counts - d1_gene_counts

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


def divideUniqueMolecules(uniqueMolecules, randomState, chromosome_division_results,
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

	# Get indexes of domains assigned to daughter 1
	d1_mother_domain_indexes = chromosome_division_results['d1_mother_domain_indexes']
	d1_all_domain_indexes = chromosome_division_results['d1_all_domain_indexes']

	# List of unique molecules that should not be binomially divided
	# Note: the only unique molecule currently not in this list is active RNA
	# polymerase.
	nonbinomial_unique_molecules = ['originOfReplication', 'chromosome_domain',
		'active_replisome', 'fullChromosome', 'activeRibosome']

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

	# Divide active replisomes according to the indexes of the chromosome
	# domains they are bound to
	moleculeSet = uniqueMolecules.container.objectsInCollection('active_replisome')
	moleculeAttributeDict = uniqueMoleculesToDivide['active_replisome']
	n_replisomes = len(moleculeSet)

	if n_replisomes > 0:
		domain_index = moleculeSet.attr("domain_index")

		# Divide replisomes based on their domain indexes
		d1_bool = np.zeros(n_replisomes, dtype=bool)
		for index in d1_all_domain_indexes:
			d1_bool = np.logical_or(d1_bool, domain_index == index)
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

		d1_unique_molecules_container.objectsNew('active_replisome', n_d1,
			**d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('active_replisome', n_d2,
			**d2_dividedAttributesDict)

	# Divide oriCs according to the indexes of the domains they belong to
	moleculeSet = uniqueMolecules.container.objectsInCollection('originOfReplication')
	moleculeAttributeDict = uniqueMoleculesToDivide['originOfReplication']
	n_oric = len(moleculeSet)

	if n_oric > 0:
		domain_index = moleculeSet.attr('domain_index')

		# Divide oriC's based on their domain indexes
		d1_bool = np.zeros(n_oric, dtype=bool)
		for index in d1_all_domain_indexes:
			d1_bool = np.logical_or(d1_bool, domain_index == index)
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

		d1_unique_molecules_container.objectsNew('originOfReplication', n_d1,
			**d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('originOfReplication', n_d2,
			**d2_dividedAttributesDict)

	# Divide chromosome domains according to their indexes
	moleculeSet = uniqueMolecules.container.objectsInCollection('chromosome_domain')
	moleculeAttributeDict = uniqueMoleculesToDivide['chromosome_domain']
	n_domain = len(moleculeSet)

	if n_domain > 0:
		domain_index = moleculeSet.attr('domain_index')

		# Divide oriC's based on their domain indexes
		d1_bool = np.zeros(n_domain, dtype=bool)
		for index in d1_all_domain_indexes:
			d1_bool = np.logical_or(d1_bool, domain_index == index)
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

		d1_unique_molecules_container.objectsNew('chromosome_domain', n_d1,
			**d1_dividedAttributesDict)
		d2_unique_molecules_container.objectsNew('chromosome_domain', n_d2,
			**d2_dividedAttributesDict)

	# Divide full chromosomes based on their mother domain indexes
	moleculeSet = uniqueMolecules.container.objectsInCollection('fullChromosome')
	moleculeAttributeDict = uniqueMoleculesToDivide['fullChromosome']
	n_full_chromosome = len(moleculeSet)

	if n_full_chromosome > 0:
		mother_domain_index = moleculeSet.attr('mother_domain_index')

		# Divide full chromosomes based on their mother domain index
		d1_bool = np.zeros(n_full_chromosome, dtype=bool)
		for index in d1_mother_domain_indexes:
			d1_bool = np.logical_or(d1_bool, mother_domain_index == index)
		d2_bool = np.logical_not(d1_bool)

		# Add the divided full chromosomes to the daughter cell containers
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


def flatten(l):
	"""
	Flattens a nested list into a single list.
	"""
	return [item for sublist in l for item in sublist]


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