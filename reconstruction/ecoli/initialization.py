"""
Functions to initialize molecule containers from sim_data.
"""

from typing import cast

import numpy as np
import scipy.sparse

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.containers.unique_objects_container import Access
from wholecell.utils.fitting import normalize, countsFromMassAndExpression, masses_and_counts_for_homeostatic_target
from wholecell.utils.polymerize import computeMassIncrease
from wholecell.utils import units
from wholecell.utils.mc_complexation import mccFormComplexesWithPrebuiltMatrices
from wholecell.utils.random import stochasticRound


RAND_MAX = 2**31


def create_bulk_container(sim_data, n_seeds=1, condition=None, seed=0, ppgpp_regulation=True, trna_attenuation=True, mass_coeff=1):
	try:
		old_condition = sim_data.condition
		if condition is not None:
			sim_data.condition = condition
		average_container = BulkObjectsContainer(sim_data.internal_state.bulk_molecules.bulk_data['id'], dtype=float)

		media_id = sim_data.conditions[sim_data.condition]['nutrients']
		exchange_data = sim_data.external_state.exchange_data_from_media(media_id)
		import_molecules = (set(exchange_data['importUnconstrainedExchangeMolecules'])
			| set(exchange_data['importConstrainedExchangeMolecules']))

		for n in range(n_seeds):
			cnt = BulkObjectsContainer(sim_data.internal_state.bulk_molecules.bulk_data['id'])
			random_state = np.random.RandomState(seed=seed+n)
			initializeBulkMolecules(cnt, sim_data, media_id, import_molecules, random_state,
				mass_coeff, ppgpp_regulation, trna_attenuation)
			average_container.countsIs(average_container.counts() + cnt.counts())
	except Exception:
		raise RuntimeError('sim_data might not be fully initialized.  Make sure all attributes have been set before using this function.')

	sim_data.condition = old_condition
	average_container.countsIs(average_container.counts() / n_seeds)
	return average_container

def initializeBulkMolecules(bulkMolCntr, sim_data, media_id, import_molecules, randomState, massCoeff,
		ppgpp_regulation, trna_attenuation):

	# Set protein counts from expression
	initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff, ppgpp_regulation, trna_attenuation)

	# Set RNA counts from expression
	initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff, ppgpp_regulation, trna_attenuation)

	# Set other biomass components
	set_small_molecule_counts(bulkMolCntr, sim_data, media_id, import_molecules, massCoeff)

	# Form complexes
	initializeComplexation(bulkMolCntr, sim_data, randomState)

def initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState,
		superhelical_density, ppgpp_regulation, trna_attenuation):
	# Initialize counts of full chromosomes
	initializeFullChromosome(bulkMolCntr, uniqueMolCntr, sim_data)

	# Initialize unique molecules relevant to replication
	initializeReplication(bulkMolCntr, uniqueMolCntr, sim_data)

	# Initialize bound transcription factors
	initialize_transcription_factors(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

	# Initialize active RNAPs and unique molecule representations of RNAs
	initialize_transcription(bulkMolCntr, uniqueMolCntr, sim_data, randomState, ppgpp_regulation, trna_attenuation)

	# Initialize linking numbers of chromosomal segments
	if superhelical_density:
		initialize_chromosomal_segments(uniqueMolCntr, sim_data)

	# Initialize activate ribosomes
	initialize_translation(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

def calculate_cell_mass(states):
	"""
	Determines the total cell mass from the currently initialized states.

	Args:
		states (dict with internal_state objects as values): internal states of sim

	Returns:
		float with mass units: total mass of the current cell state
	"""

	mass = 0
	for state in states.values():
		state.calculateMass()
		mass += np.sum(state.mass())

	return units.fg * mass

def initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff, ppgpp_regulation, trna_attenuation):

	monomersView = bulkMolCntr.countsView(sim_data.process.translation.monomer_data["id"])
	monomerMass = massCoeff * sim_data.mass.get_component_masses(sim_data.condition_to_doubling_time[sim_data.condition])["proteinMass"] / sim_data.mass.avg_cell_to_initial_cell_conversion_factor
	# TODO: unify this logic with the parca so it doesn't fall out of step
	# again (look at the calcProteinCounts function)

	transcription = sim_data.process.transcription
	if ppgpp_regulation:
		rnaExpression = sim_data.calculate_ppgpp_expression(sim_data.condition)
	else:
		rnaExpression = transcription.rna_expression[sim_data.condition]

	if trna_attenuation:
		# Need to adjust expression (calculated without attenuation) by basal_adjustment
		# to get the expected expression without any attenuation and then multiply
		# by the condition readthrough probability to get the condition specific expression
		readthrough = transcription.attenuation_readthrough[sim_data.condition]
		basal_adjustment = transcription.attenuation_readthrough['basal']
		rnaExpression[transcription.attenuated_rna_indices] *= readthrough / basal_adjustment

	monomerExpression = normalize(
		sim_data.process.transcription.cistron_tu_mapping_matrix.dot(rnaExpression)[
			sim_data.relation.cistron_to_monomer_mapping] *
		sim_data.process.translation.translation_efficiencies_by_monomer /
		(np.log(2) / sim_data.condition_to_doubling_time[sim_data.condition].asNumber(units.s) + sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s))
		)

	nMonomers = countsFromMassAndExpression(
		monomerMass.asNumber(units.g),
		sim_data.process.translation.monomer_data["mw"].asNumber(units.g / units.mol),
		monomerExpression,
		sim_data.constants.n_avogadro.asNumber(1 / units.mol)
		)

	monomersView.countsIs(
		randomState.multinomial(nMonomers, monomerExpression)
		)


def initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff, ppgpp_regulation, trna_attenuation):
	"""
	Initializes counts of RNAs in the bulk molecule container using RNA
	expression data. mRNA counts are also initialized here, but is later reset
	to zero when the representations for mRNAs are moved to the unique molecule
	container.
	"""

	transcription = sim_data.process.transcription

	rnaView = bulkMolCntr.countsView(transcription.rna_data["id"])
	rnaMass = massCoeff * sim_data.mass.get_component_masses(sim_data.condition_to_doubling_time[sim_data.condition])["rnaMass"] / sim_data.mass.avg_cell_to_initial_cell_conversion_factor

	if ppgpp_regulation:
		rnaExpression = sim_data.calculate_ppgpp_expression(sim_data.condition)
	else:
		rnaExpression = normalize(transcription.rna_expression[sim_data.condition])

	if trna_attenuation:
		# Need to adjust expression (calculated without attenuation) by basal_adjustment
		# to get the expected expression without any attenuation and then multiply
		# by the condition readthrough probability to get the condition specific expression
		readthrough = transcription.attenuation_readthrough[sim_data.condition]
		basal_adjustment = transcription.attenuation_readthrough['basal']
		rnaExpression[transcription.attenuated_rna_indices] *= readthrough / basal_adjustment
		rnaExpression /= rnaExpression.sum()

	nRnas = countsFromMassAndExpression(
		rnaMass.asNumber(units.g),
		transcription.rna_data["mw"].asNumber(units.g / units.mol),
		rnaExpression,
		sim_data.constants.n_avogadro.asNumber(1 / units.mol)
		)

	# ID Groups of rRNAs
	idx_16Srrna = np.where(transcription.rna_data['is_16S_rRNA'])[0]
	idx_23Srrna = np.where(transcription.rna_data['is_23S_rRNA'])[0]
	idx_5Srrna = np.where(transcription.rna_data['is_5S_rRNA'])[0]

	# Assume expression from all rRNA genes produce rRNAs from the first operon
	total_16Srrna_expression = rnaExpression[idx_16Srrna].sum()
	total_23Srrna_expression = rnaExpression[idx_23Srrna].sum()
	total_5Srrna_expression = rnaExpression[idx_5Srrna].sum()

	rnaExpression[idx_16Srrna] = 0
	rnaExpression[idx_23Srrna] = 0
	rnaExpression[idx_5Srrna] = 0

	rnaExpression[idx_16Srrna[0]] = total_16Srrna_expression
	rnaExpression[idx_23Srrna[0]] = total_23Srrna_expression
	rnaExpression[idx_5Srrna[0]] = total_5Srrna_expression

	# Calculate initial counts of each RNA from multinomial distribution
	rnaView.countsIs(
		randomState.multinomial(nRnas, rnaExpression)
		)

# TODO: remove checks for zero concentrations (change to assertion)
# TODO: move any rescaling logic to KB/fitting
def set_small_molecule_counts(bulkMolCntr, sim_data, media_id, import_molecules, massCoeff, cell_mass=None):
	doubling_time = sim_data.condition_to_doubling_time[sim_data.condition]

	concDict = sim_data.process.metabolism.concentration_updates.concentrations_based_on_nutrients(
		media_id=media_id, imports=import_molecules
		)
	concDict.update(sim_data.mass.getBiomassAsConcentrations(doubling_time))
	concDict[sim_data.molecule_ids.ppGpp] = sim_data.growth_rate_parameters.get_ppGpp_conc(doubling_time)
	moleculeIds = sorted(concDict)
	moleculeConcentrations = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in moleculeIds])

	if cell_mass is None:
		avgCellFractionMass = sim_data.mass.get_component_masses(doubling_time)
		other_dry_mass = (massCoeff * (avgCellFractionMass["proteinMass"]
			+ avgCellFractionMass["rnaMass"] + avgCellFractionMass["dnaMass"])
			/ sim_data.mass.avg_cell_to_initial_cell_conversion_factor)
	else:
		small_molecule_mass = 0 * units.fg
		for mol in concDict:
			small_molecule_mass += (bulkMolCntr.count(mol)
									* sim_data.getter.get_mass(mol) / sim_data.constants.n_avogadro)
		other_dry_mass = cell_mass - small_molecule_mass

	massesToAdd, countsToAdd = masses_and_counts_for_homeostatic_target(
		other_dry_mass,
		moleculeConcentrations,
		sim_data.getter.get_masses(moleculeIds),
		sim_data.constants.cell_density,
		sim_data.constants.n_avogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		moleculeIds
		)

def initializeComplexation(bulkMolCntr, sim_data, randomState):
	moleculeNames = sim_data.process.complexation.molecule_names
	moleculeView = bulkMolCntr.countsView(moleculeNames)

	# save rnase counts for after complexation, need to be monomers
	rnases = sim_data.process.rna_decay.endoRNase_ids
	rnaseCounts = bulkMolCntr.countsView(rnases).counts()
	bulkMolCntr.countsIs(0, rnases)

	stoichMatrix = sim_data.process.complexation.stoich_matrix().astype(np.int64, order='F')

	moleculeCounts = moleculeView.counts()
	updatedMoleculeCounts, complexationEvents = mccFormComplexesWithPrebuiltMatrices(
		moleculeCounts,
		randomState.randint(1000),
		stoichMatrix,
		*sim_data.process.complexation.prebuilt_matrices)

	bulkMolCntr.countsIs(
		updatedMoleculeCounts,
		moleculeNames)

	if np.any(updatedMoleculeCounts < 0):
		raise ValueError('Negative counts after complexation')

	bulkMolCntr.countsIs(rnaseCounts, rnases)


def initializeFullChromosome(bulkMolCntr, uniqueMolCntr, sim_data):
	"""
	Initializes the counts of full chromosomes to one. The division_time of
	this initial chromosome is set to be zero for consistency.
	"""
	uniqueMolCntr.objectsNew(
		'full_chromosome', 1,
		division_time=0.0,
		has_triggered_division=True,
		domain_index=0,
		)


def initializeReplication(bulkMolCntr, uniqueMolCntr, sim_data):
	"""
	Initializes replication by creating an appropriate number of replication
	forks given the cell growth rate. This also initializes the gene dosage
	bulk counts using the initial locations of the forks.
	"""
	# Determine the number and location of replication forks at the start of
	# the cell cycle
	# Get growth rate constants
	C = sim_data.process.replication.c_period
	D = sim_data.process.replication.d_period
	tau = sim_data.condition_to_doubling_time[sim_data.condition].asUnit(units.min)

	# Calculate length of replichore
	genome_length = sim_data.process.replication.genome_length
	replichore_length = np.ceil(0.5*genome_length) * units.nt

	# Calculate the maximum number of replisomes that could be formed with
	# the existing counts of replisome subunits
	n_max_replisomes = np.min(np.concatenate(
		(bulkMolCntr.counts(sim_data.molecule_groups.replisome_trimer_subunits) // 3,
		bulkMolCntr.counts(sim_data.molecule_groups.replisome_monomer_subunits))))

	# Generate arrays specifying appropriate initial replication conditions
	oric_state, replisome_state, domain_state = determine_chromosome_state(
		C, D, tau, replichore_length, n_max_replisomes,
		sim_data.process.replication.no_child_place_holder)

	n_oric = oric_state["domain_index"].size
	n_replisome = replisome_state["domain_index"].size
	n_domain = domain_state["domain_index"].size

	# Add OriC molecules with the proposed attributes
	uniqueMolCntr.objectsNew(
		'oriC', n_oric,
		domain_index=oric_state["domain_index"])

	# Add chromosome domain molecules with the proposed attributes
	uniqueMolCntr.objectsNew(
		"chromosome_domain", n_domain,
		domain_index=domain_state["domain_index"],
		child_domains=domain_state["child_domains"])

	if n_replisome != 0:
		# Update mass to account for DNA strands that have already been
		# elongated.
		sequences = sim_data.process.replication.replication_sequences
		fork_coordinates = replisome_state["coordinates"]
		sequence_elongations = np.abs(np.repeat(fork_coordinates, 2))

		mass_increase_dna = computeMassIncrease(
			np.tile(sequences, (n_replisome//2, 1)),
			sequence_elongations,
			sim_data.process.replication.replication_monomer_weights.asNumber(units.fg))

		# Add active replisomes as unique molecules and set attributes
		uniqueMolCntr.objectsNew(
			'active_replisome', n_replisome,
			coordinates=replisome_state["coordinates"],
			right_replichore=replisome_state["right_replichore"],
			domain_index=replisome_state["domain_index"],
			massDiff_DNA=mass_increase_dna[0::2] + mass_increase_dna[1::2])

		# Remove replisome subunits from bulk molecules
		bulkMolCntr.countsDec(3*n_replisome, sim_data.molecule_groups.replisome_trimer_subunits)
		bulkMolCntr.countsDec(n_replisome, sim_data.molecule_groups.replisome_monomer_subunits)

	# Get coordinates of all promoters and DnaA boxes
	all_promoter_coordinates = sim_data.process.transcription.rna_data[
		"replication_coordinate"]
	all_DnaA_box_coordinates = sim_data.process.replication.motif_coordinates[
		"DnaA_box"]


	# Define function that initializes attributes of sequence motifs given the
	# initial state of the chromosome
	def get_motif_attributes(all_motif_coordinates):
		"""
		Using the initial positions of replication forks, calculate attributes
		of unique molecules representing DNA motifs, given their positions on
		the genome.

		Args:
			all_motif_coordinates (ndarray): Genomic coordinates of DNA motifs,
			represented in a specific order

		Returns:
			motif_index: Indices of all motif copies, in the case where
			different indexes imply a different functional role
			motif_coordinates: Genomic coordinates of all motif copies
			motif_domain_index: Domain indexes of the chromosome domain that
			each motif copy belongs to
		"""
		motif_index, motif_coordinates, motif_domain_index = [], [], []

		def in_bounds(coordinates, lb, ub):
			return np.logical_and(coordinates < ub, coordinates > lb)

		# Loop through all chromosome domains
		for domain_idx in domain_state["domain_index"]:

			# If the domain is the mother domain of the initial chromosome,
			if domain_idx == 0:
				if n_replisome == 0:
					# No replisomes - all motifs should fall in this domain
					motif_mask = np.ones_like(all_motif_coordinates, dtype=bool)

				else:
					# Get domain boundaries
					domain_boundaries = replisome_state["coordinates"][
						replisome_state["domain_index"] == 0]

					# Add motifs outside of this boundary
					motif_mask = np.logical_or(
						all_motif_coordinates > domain_boundaries.max(),
						all_motif_coordinates < domain_boundaries.min())

			# If the domain contains the origin,
			elif np.isin(domain_idx, oric_state["domain_index"]):
				# Get index of the parent domain
				parent_domain_idx = domain_state["domain_index"][
					np.where(domain_state["child_domains"] == domain_idx)[0]]

				# Get domain boundaries of the parent domain
				parent_domain_boundaries = replisome_state["coordinates"][
					replisome_state["domain_index"] == parent_domain_idx]

				# Add motifs inside this boundary
				motif_mask = in_bounds(all_motif_coordinates,
					parent_domain_boundaries.min(),
					parent_domain_boundaries.max())

			# If the domain neither contains the origin nor the terminus,
			else:
				# Get index of the parent domain
				parent_domain_idx = domain_state["domain_index"][
					np.where(domain_state["child_domains"] == domain_idx)[0]]

				# Get domain boundaries of the parent domain
				parent_domain_boundaries = replisome_state["coordinates"][
					replisome_state["domain_index"] == parent_domain_idx]

				# Get domain boundaries of this domain
				domain_boundaries = replisome_state["coordinates"][
					replisome_state["domain_index"] == domain_idx]

				# Add motifs between the boundaries
				motif_mask = np.logical_or(
					in_bounds(all_motif_coordinates,
						domain_boundaries.max(),
						parent_domain_boundaries.max()),
					in_bounds(all_motif_coordinates,
						parent_domain_boundaries.min(),
						domain_boundaries.min())
					)

			# Append attributes to existing list
			motif_index.extend(np.nonzero(motif_mask)[0])
			motif_coordinates.extend(all_motif_coordinates[motif_mask])
			motif_domain_index.extend(np.full(motif_mask.sum(), domain_idx))

		return motif_index, motif_coordinates, motif_domain_index


	# Use function to get attributes for promoters and DnaA boxes
	TU_index, promoter_coordinates, promoter_domain_index = get_motif_attributes(
		all_promoter_coordinates)
	_, DnaA_box_coordinates, DnaA_box_domain_index = get_motif_attributes(
		all_DnaA_box_coordinates)

	# Add promoters as unique molecules and set attributes
	# Note: the bound_TF attribute is properly initialized in the function
	# initialize_transcription_factors
	n_promoter = len(TU_index)
	n_tf = len(sim_data.process.transcription_regulation.tf_ids)

	uniqueMolCntr.objectsNew(
		'promoter', n_promoter,
		TU_index=np.array(TU_index),
		coordinates=np.array(promoter_coordinates),
		domain_index=np.array(promoter_domain_index),
		bound_TF=np.zeros((n_promoter, n_tf), dtype=bool))

	# Add DnaA boxes as unique molecules and set attributes
	n_DnaA_box = len(DnaA_box_coordinates)

	uniqueMolCntr.objectsNew(
		'DnaA_box', n_DnaA_box,
		coordinates=np.array(DnaA_box_coordinates),
		domain_index=np.array(DnaA_box_domain_index),
		DnaA_bound=np.zeros(n_DnaA_box, dtype=bool)
		)


def initialize_transcription_factors(bulkMolCntr, uniqueMolCntr, sim_data, randomState):
	"""
	Initialize transcription factors that are bound to the chromosome. For each
	type of transcription factor, this function calculates the total number of
	transcription factors that should be bound to the chromosome using the
	binding probabilities of each transcription factor and the number of
	available promoter sites. The calculated number of transcription factors
	are then distributed randomly to promoters, whose bound_TF attributes and
	submasses are updated correspondingly.
	"""
	# Get transcription factor properties from sim_data
	tf_ids = sim_data.process.transcription_regulation.tf_ids
	tfToTfType = sim_data.process.transcription_regulation.tf_to_tf_type
	pPromoterBoundTF = sim_data.process.transcription_regulation.p_promoter_bound_tf

	# Build dict that maps TFs to transcription units they regulate
	delta_prob = sim_data.process.transcription_regulation.delta_prob
	TF_to_TU_idx = {}

	for i, tf in enumerate(tf_ids):
		TF_to_TU_idx[tf] = delta_prob['deltaI'][delta_prob['deltaJ'] == i]

	# Get views into bulk molecule representations of transcription factors
	active_tf_view = {}
	inactive_tf_view = {}

	for tf in tf_ids:
		active_tf_view[tf] = bulkMolCntr.countsView([tf + "[c]"])

		if tfToTfType[tf] == "1CS":
			if tf == sim_data.process.transcription_regulation.active_to_bound[tf]:
				inactive_tf_view[tf] = bulkMolCntr.countsView([
					sim_data.process.equilibrium.get_unbound(tf + "[c]")
					])
			else:
				inactive_tf_view[tf] = bulkMolCntr.countsView([
					sim_data.process.transcription_regulation.active_to_bound[tf] + "[c]"
					])
		elif tfToTfType[tf] == "2CS":
			inactive_tf_view[tf] = bulkMolCntr.countsView([
				sim_data.process.two_component_system.active_to_inactive_tf[tf + "[c]"]
				])

	# Get masses of active transcription factors
	bulk_molecule_ids = sim_data.internal_state.bulk_molecules.bulk_data["id"]
	tf_indexes = [np.where(bulk_molecule_ids == tf_id + "[c]")[0][0]
		for tf_id in tf_ids]
	active_tf_masses = (sim_data.internal_state.bulk_molecules.bulk_data["mass"][
		tf_indexes] / sim_data.constants.n_avogadro).asNumber(units.fg)

	# Get attributes of promoters
	promoters = uniqueMolCntr.objectsInCollection(
		'promoter', access=[Access.EDIT])
	TU_index = promoters.attr("TU_index")

	# Initialize bound_TF array
	bound_TF = np.zeros((len(promoters), len(tf_ids)), dtype=bool)

	for tf_idx, tf_id in enumerate(tf_ids):
		# Get counts of transcription factors
		active_tf_counts = active_tf_view[tf_id].counts()

		# If there are no active transcription factors at initialization,
		# continue to the next transcription factor
		if active_tf_counts == 0:
			continue

		# Compute probability of binding the promoter
		if tfToTfType[tf_id] == "0CS":
			pPromoterBound = 1.
		else:
			inactive_tf_counts = inactive_tf_view[tf_id].counts()
			pPromoterBound = pPromoterBoundTF(
				active_tf_counts, inactive_tf_counts)

		# Determine the number of available promoter sites
		available_promoters = np.isin(TU_index, TF_to_TU_idx[tf_id])
		n_available_promoters = available_promoters.sum()

		# Calculate the number of promoters that should be bound
		n_to_bind = int(stochasticRound(
			randomState, np.full(n_available_promoters, pPromoterBound)).sum())

		bound_locs = np.zeros(n_available_promoters, dtype=bool)
		if n_to_bind > 0:
			# Determine randomly which DNA targets to bind based on which of
			# the following is more limiting:
			# number of promoter sites to bind, or number of active
			# transcription factors
			bound_locs[
				randomState.choice(
					n_available_promoters,
					size=min(n_to_bind, active_tf_view[tf_id].counts()[0]),
					replace=False)
			] = True

			# Update count of free transcription factors
			active_tf_view[tf_id].countsDec(bound_locs.sum())

			# Update bound_TF array
			bound_TF[available_promoters, tf_idx] = bound_locs

	# Calculate masses of bound TFs
	mass_diffs = bound_TF.dot(active_tf_masses)

	# Reset bound_TF attribute of promoters
	promoters.attrIs(bound_TF=bound_TF)

	# Add mass_diffs array to promoter submass
	promoters.add_submass_by_array(mass_diffs)


def initialize_transcription(bulkMolCntr, uniqueMolCntr, sim_data, randomState,
		ppgpp_regulation, trna_attenuation):
	"""
	Activate RNA polymerases as unique molecules, and distribute them along
	lengths of trancription units, while decreasing counts of unactivated RNA
	polymerases (APORNAP-CPLX[c]). Also initialize unique molecule
	representations of fully transcribed mRNAs and partially transcribed RNAs,
	using counts of mRNAs initialized as bulk molecules, and the attributes of
	initialized RNA polymerases. The counts of full mRNAs represented as bulk
	molecules are reset to zero.

	RNA polymerases are placed randomly across the length of each transcription
	unit, with the synthesis probabilities for each TU determining the number of
	RNA polymerases placed at each gene.
	"""
	# Load parameters
	rnaLengths = sim_data.process.transcription.rna_data['length'].asNumber()
	rna_masses = (sim_data.process.transcription.rna_data['mw'] / sim_data.constants.n_avogadro).asNumber(units.fg)
	current_media_id = sim_data.conditions[sim_data.condition]['nutrients']
	fracActiveRnap = sim_data.process.transcription.rnapFractionActiveDict[current_media_id]
	inactive_RNAP_counts = bulkMolCntr.countsView([sim_data.molecule_ids.full_RNAP]).counts()[0]
	rnaSequences = sim_data.process.transcription.transcription_sequences
	ntWeights = sim_data.process.transcription.transcription_monomer_weights
	endWeight = sim_data.process.transcription.transcription_end_weight
	replichore_lengths = sim_data.process.replication.replichore_lengths
	chromosome_length = replichore_lengths.sum()

	# Number of rnaPoly to activate
	n_RNAPs_to_activate = np.int64(fracActiveRnap * inactive_RNAP_counts)

	# Get attributes of promoters
	promoters = uniqueMolCntr.objectsInCollection("promoter")
	n_promoters = len(promoters)
	TU_index, bound_TF, domain_index_promoters = promoters.attrs(
		"TU_index", "bound_TF", "domain_index")

	# Parameters for rnaSynthProb
	if ppgpp_regulation:
		doubling_time = sim_data.condition_to_doubling_time[sim_data.condition]
		ppgpp_conc = sim_data.growth_rate_parameters.get_ppGpp_conc(doubling_time)
		basal_prob, _= sim_data.process.transcription.synth_prob_from_ppgpp(
			ppgpp_conc, sim_data.process.replication.get_average_copy_number)
		ppgpp_scale = basal_prob[TU_index]
		ppgpp_scale[ppgpp_scale == 0] = 1  # Use original delta prob if no ppGpp basal prob
	else:
		basal_prob = sim_data.process.transcription_regulation.basal_prob.copy()
		ppgpp_scale = 1

	if trna_attenuation:
		basal_prob[sim_data.process.transcription.attenuated_rna_indices] += sim_data.process.transcription.attenuation_basal_prob_adjustments
	n_TUs = len(basal_prob)
	delta_prob_matrix = sim_data.process.transcription_regulation.get_delta_prob_matrix(dense=True, ppgpp=ppgpp_regulation)

	# Construct matrix that maps promoters to transcription units
	TU_to_promoter = scipy.sparse.csr_matrix(
		(np.ones(n_promoters), (TU_index, np.arange(n_promoters))),
		shape=(n_TUs, n_promoters))

	# Synthesis probabilities for different categories of genes
	rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction
	rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
	rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase

	# Get coordinates and transcription directions of transcription units
	replication_coordinate = sim_data.process.transcription.rna_data[
		"replication_coordinate"]
	transcription_direction = sim_data.process.transcription.rna_data[
		"direction"]

	# Determine changes from genetic perturbations
	genetic_perturbations = {}
	perturbations = getattr(sim_data, 'genetic_perturbations', {})

	if len(perturbations) > 0:
		probability_indexes = [
			(index, sim_data.genetic_perturbations[rna_data['id']])
				for index, rna_data in enumerate(sim_data.process.transcription.rna_data)
				if rna_data['id'] in sim_data.genetic_perturbations]

		genetic_perturbations = {
			'fixedRnaIdxs': [pair[0] for pair in probability_indexes],
			'fixedSynthProbs': [pair[1] for pair in probability_indexes]}

	# If initiationShuffleIdxs does not exist, set value to None
	shuffleIdxs = getattr(sim_data.process.transcription, 'initiationShuffleIdxs', None)

	# ID Groups
	idx_rRNA = np.where(sim_data.process.transcription.rna_data['is_rRNA'])[0]
	idx_mRNA = np.where(sim_data.process.transcription.rna_data['is_mRNA'])[0]
	idx_tRNA = np.where(sim_data.process.transcription.rna_data['is_tRNA'])[0]
	idx_rprotein = np.where(sim_data.process.transcription.rna_data['includes_ribosomal_protein'])[0]
	idx_rnap = np.where(sim_data.process.transcription.rna_data['includes_RNAP'])[0]

	# Calculate probabilities of the RNAP binding to the promoters
	promoter_init_probs = (basal_prob[TU_index] + ppgpp_scale *
		np.multiply(delta_prob_matrix[TU_index, :], bound_TF).sum(axis=1))

	if len(genetic_perturbations) > 0:
		rescale_initiation_probs(
			promoter_init_probs, TU_index,
			genetic_perturbations["fixedSynthProbs"],
			genetic_perturbations["fixedRnaIdxs"])

	# Adjust probabilities to not be negative
	promoter_init_probs[promoter_init_probs < 0] = 0.0
	promoter_init_probs /= promoter_init_probs.sum()
	if np.any(promoter_init_probs < 0):
		raise Exception("Have negative RNA synthesis probabilities")

	# Adjust synthesis probabilities depending on environment
	synthProbFractions = rnaSynthProbFractions[current_media_id]

	# Create masks for different types of RNAs
	is_mRNA = np.isin(TU_index, idx_mRNA)
	is_tRNA = np.isin(TU_index, idx_tRNA)
	is_rRNA = np.isin(TU_index, idx_rRNA)
	is_rprotein = np.isin(TU_index, idx_rprotein)
	is_rnap = np.isin(TU_index, idx_rnap)
	is_fixed = is_tRNA | is_rRNA | is_rprotein | is_rnap

	# Rescale initiation probabilities based on type of RNA
	promoter_init_probs[is_mRNA] *= synthProbFractions["mRna"] / promoter_init_probs[is_mRNA].sum()
	promoter_init_probs[is_tRNA] *= synthProbFractions["tRna"] / promoter_init_probs[is_tRNA].sum()
	promoter_init_probs[is_rRNA] *= synthProbFractions["rRna"] / promoter_init_probs[is_rRNA].sum()

	# Set fixed synthesis probabilities for RProteins and RNAPs
	rescale_initiation_probs(
		promoter_init_probs, TU_index,
		np.concatenate((
			rnaSynthProbRProtein[current_media_id],
			rnaSynthProbRnaPolymerase[current_media_id])),
		np.concatenate((idx_rprotein, idx_rnap)))

	assert promoter_init_probs[is_fixed].sum() < 1.0

	# Adjust for attenuation that will stop transcription after initiation
	if trna_attenuation:
		attenuation_readthrough = {
			idx: prob for idx, prob in
			zip(sim_data.process.transcription.attenuated_rna_indices, sim_data.process.transcription.attenuation_readthrough[sim_data.condition])
			}
		readthrough_adjustment = np.array([attenuation_readthrough.get(idx, 1) for idx in TU_index])
		promoter_init_probs *= readthrough_adjustment

	scaleTheRestBy = (1. - promoter_init_probs[is_fixed].sum()) / promoter_init_probs[~is_fixed].sum()
	promoter_init_probs[~is_fixed] *= scaleTheRestBy

	# Compute synthesis probabilities of each transcription unit
	TU_synth_probs = TU_to_promoter.dot(promoter_init_probs)

	# Shuffle initiation rates if we're running the variant that calls this
	if shuffleIdxs is not None:
		rescale_initiation_probs(
			promoter_init_probs, TU_index, TU_synth_probs[shuffleIdxs],
			np.arange(n_TUs))

	# normalize to length of rna
	init_prob_length_adjusted = promoter_init_probs * rnaLengths[TU_index]
	init_prob_normalized = init_prob_length_adjusted / init_prob_length_adjusted.sum()

	# Sample a multinomial distribution of synthesis probabilities to determine
	# what RNA are initialized
	n_initiations = randomState.multinomial(
		n_RNAPs_to_activate, init_prob_normalized)

	# Build array of transcription unit indexes for partially transcribed mRNAs
	# and domain indexes for RNAPs
	TU_index_partial_RNAs = np.repeat(TU_index, n_initiations)
	domain_index_rnap = np.repeat(domain_index_promoters, n_initiations)

	# Build arrays of starting coordinates and transcription directions
	starting_coordinates = replication_coordinate[TU_index_partial_RNAs]
	direction = transcription_direction[TU_index_partial_RNAs]

	# Randomly advance RNAPs along the transcription units
	# TODO (Eran): make sure there aren't any RNAPs at same location on same TU
	updated_lengths = np.array(
		randomState.rand(n_RNAPs_to_activate) * rnaLengths[TU_index_partial_RNAs],
		dtype=int)

	# Rescale boolean array of directions to an array of 1's and -1's.
	direction_rescaled = (2 * (direction - 0.5)).astype(np.int64)

	# Compute the updated coordinates of RNAPs. Coordinates of RNAPs moving in
	# the positive direction are increased, whereas coordinates of RNAPs moving
	# in the negative direction are decreased.
	updated_coordinates = starting_coordinates + np.multiply(
		direction_rescaled, updated_lengths)

	# Reset coordinates of RNAPs that cross the boundaries between right and
	# left replichores
	updated_coordinates[
		updated_coordinates > replichore_lengths[0]] -= chromosome_length
	updated_coordinates[
		updated_coordinates < -replichore_lengths[1]] += chromosome_length

	# Update mass
	sequences = rnaSequences[TU_index_partial_RNAs]
	added_mass = computeMassIncrease(sequences, updated_lengths, ntWeights)
	added_mass[updated_lengths != 0] += endWeight  # add endWeight to all new Rna

	# Masses of partial mRNAs are counted as mRNA mass as they are already
	# functional, but the masses of other types of partial RNAs are counted as
	# generic RNA mass.
	added_RNA_mass = added_mass.copy()
	added_mRNA_mass = added_mass.copy()

	is_mRNA_partial_RNAs = np.isin(TU_index_partial_RNAs, idx_mRNA)
	added_RNA_mass[is_mRNA_partial_RNAs] = 0
	added_mRNA_mass[np.logical_not(is_mRNA_partial_RNAs)] = 0

	# Add active RNAPs and get their unique indexes
	RNAP_indexes = uniqueMolCntr.objectsNew(
		'active_RNAP', n_RNAPs_to_activate,
		domain_index=domain_index_rnap,
		coordinates=updated_coordinates,
		direction=direction)

	# Decrement counts of bulk inactive RNAPs
	bulkMolCntr.countsIs(inactive_RNAP_counts - n_RNAPs_to_activate,
		[sim_data.molecule_ids.full_RNAP])

	# Add partially transcribed RNAs
	uniqueMolCntr.objectsNew(
		'RNA', n_RNAPs_to_activate,
		TU_index=TU_index_partial_RNAs,
		transcript_length=updated_lengths,
		is_mRNA=is_mRNA_partial_RNAs,
		is_full_transcript=np.zeros(cast(int, n_RNAPs_to_activate), dtype=bool),
		can_translate=is_mRNA_partial_RNAs,
		RNAP_index=RNAP_indexes,
		massDiff_nonspecific_RNA=added_RNA_mass,
		massDiff_mRNA=added_mRNA_mass)

	# Get counts of mRNAs initialized as bulk molecules
	mRNA_ids = sim_data.process.transcription.rna_data["id"][
		sim_data.process.transcription.rna_data['is_mRNA']]
	mRNA_view = bulkMolCntr.countsView(mRNA_ids)
	mRNA_counts = mRNA_view.counts()

	# Subtract number of partially transcribed mRNAs that were initialized.
	# Note: some mRNAs with high degradation rates have more partial mRNAs than
	# the expected total number of mRNAs - for these mRNAs we simply set the
	# initial full mRNA counts to be zero.
	partial_mRNA_counts = np.bincount(
		TU_index_partial_RNAs[is_mRNA_partial_RNAs], minlength=n_TUs)[idx_mRNA]
	full_mRNA_counts = (mRNA_counts - partial_mRNA_counts).clip(min=0)

	# Get array of TU indexes for each full mRNA
	TU_index_full_mRNAs = np.repeat(idx_mRNA, full_mRNA_counts)

	# Add fully transcribed mRNAs. The RNAP_index attribute of these molecules
	# are set to -1.
	uniqueMolCntr.objectsNew(
		'RNA', len(TU_index_full_mRNAs),
		TU_index=TU_index_full_mRNAs,
		transcript_length=rnaLengths[TU_index_full_mRNAs],
		is_mRNA=np.ones_like(TU_index_full_mRNAs, dtype=bool),
		is_full_transcript=np.ones_like(TU_index_full_mRNAs, dtype=bool),
		can_translate=np.ones_like(TU_index_full_mRNAs, dtype=bool),
		RNAP_index=np.full(TU_index_full_mRNAs.shape, -1, dtype=np.int64),
		massDiff_mRNA=rna_masses[TU_index_full_mRNAs])

	# Reset counts of bulk mRNAs to zero
	mRNA_view.countsIs(0)


def initialize_chromosomal_segments(uniqueMolCntr, sim_data):
	"""
	Initialize unique molecule representations of chromosomal segments. All
	chromosomal segments are assumed to be at their relaxed states upon
	initialization.
	"""
	# Load parameters
	relaxed_DNA_base_pairs_per_turn = sim_data.process.chromosome_structure.relaxed_DNA_base_pairs_per_turn
	terC_index = sim_data.process.chromosome_structure.terC_dummy_molecule_index
	replichore_lengths = sim_data.process.replication.replichore_lengths
	min_coordinates = -replichore_lengths[1]
	max_coordinates = replichore_lengths[0]

	# Get attributes of replisomes, active RNAPs, chromosome domains, full
	# chromosomes, and oriCs
	replisomes = uniqueMolCntr.objectsInCollection('active_replisome')
	if len(replisomes) > 0:
		replisome_coordinates, replisome_domain_indexes, replisome_unique_indexes = replisomes.attrs(
			'coordinates', 'domain_index', 'unique_index')
	else:
		replisome_coordinates = np.array([])
		replisome_domain_indexes = np.array([])
		replisome_unique_indexes = np.array([])

	active_RNAPs = uniqueMolCntr.objectsInCollection('active_RNAP')
	active_RNAP_coordinates, active_RNAP_domain_indexes, active_RNAP_unique_indexes = active_RNAPs.attrs(
		'coordinates', 'domain_index', 'unique_index')

	chromosome_domains = uniqueMolCntr.objectsInCollection('chromosome_domain')
	chromosome_domain_domain_indexes, child_domains = chromosome_domains.attrs(
		'domain_index', 'child_domains')

	full_chromosomes = uniqueMolCntr.objectsInCollection('full_chromosome')
	full_chromosome_domain_indexes = full_chromosomes.attr('domain_index')

	oriCs = uniqueMolCntr.objectsInCollection('oriC')
	origin_domain_indexes = oriCs.attr('domain_index')

	# Initialize chromosomal segment attributes
	all_boundary_molecule_indexes = np.empty((0, 2), dtype=np.int64)
	all_boundary_coordinates = np.empty((0, 2), dtype=np.int64)
	all_segment_domain_indexes = np.array([], dtype=np.int32)
	all_linking_numbers = np.array([], dtype=np.float64)


	def get_chromosomal_segment_attributes(coordinates, unique_indexes,
			spans_oriC, spans_terC):
		"""
		Returns the attributes of all chromosomal segments from a continuous
		stretch of DNA, given the coordinates and unique indexes of all
		boundary molecules.
		"""
		coordinates_argsort = np.argsort(coordinates)
		coordinates_sorted = coordinates[coordinates_argsort]
		unique_indexes_sorted = unique_indexes[coordinates_argsort]

		# Add dummy molecule at terC if domain spans terC
		if spans_terC:
			coordinates_sorted = np.insert(
				coordinates_sorted, [0, len(coordinates_sorted)],
				[min_coordinates, max_coordinates])
			unique_indexes_sorted = np.insert(
				unique_indexes_sorted, [0, len(unique_indexes_sorted)],
				terC_index)

		boundary_molecule_indexes = np.hstack((
			unique_indexes_sorted[:-1][:, np.newaxis],
			unique_indexes_sorted[1:][:, np.newaxis]))
		boundary_coordinates = np.hstack((
			coordinates_sorted[:-1][:, np.newaxis],
			coordinates_sorted[1:][:, np.newaxis]))

		# Remove segment that spans oriC if the domain does not span oriC
		if not spans_oriC:
			oriC_segment_index = np.where(
				np.sign(boundary_coordinates).sum(axis=1) == 0)[0]
			assert len(oriC_segment_index) == 1

			boundary_molecule_indexes = np.delete(boundary_molecule_indexes,
				oriC_segment_index, 0)
			boundary_coordinates = np.delete(boundary_coordinates,
				oriC_segment_index, 0)

		# Assumes all segments are at their relaxed state at initialization
		linking_numbers = (
			boundary_coordinates[:, 1] - boundary_coordinates[:, 0]
			) / relaxed_DNA_base_pairs_per_turn

		return boundary_molecule_indexes, boundary_coordinates, linking_numbers


	# Loop through each domain index
	for domain_index in chromosome_domain_domain_indexes:
		domain_spans_oriC = (domain_index in origin_domain_indexes)
		domain_spans_terC = (domain_index in full_chromosome_domain_indexes)

		# Get coordinates and indexes of all RNAPs on this domain
		RNAP_domain_mask = (active_RNAP_domain_indexes == domain_index)
		molecule_coordinates_this_domain = active_RNAP_coordinates[RNAP_domain_mask]
		molecule_indexes_this_domain = active_RNAP_unique_indexes[RNAP_domain_mask]

		# Append coordinates and indexes of replisomes on this domain, if any
		if not domain_spans_oriC:
			replisome_domain_mask = (replisome_domain_indexes == domain_index)
			molecule_coordinates_this_domain = np.concatenate((
				molecule_coordinates_this_domain,
				replisome_coordinates[replisome_domain_mask]))
			molecule_indexes_this_domain = np.concatenate((
				molecule_indexes_this_domain,
				replisome_unique_indexes[replisome_domain_mask]
				))

		# Append coordinates and indexes of parent domain replisomes, if any
		if not domain_spans_terC:
			parent_domain_index = chromosome_domain_domain_indexes[
				np.where(child_domains == domain_index)[0][0]]
			replisome_parent_domain_mask = (replisome_domain_indexes == parent_domain_index)
			molecule_coordinates_this_domain = np.concatenate((
				molecule_coordinates_this_domain,
				replisome_coordinates[replisome_parent_domain_mask]))
			molecule_indexes_this_domain = np.concatenate((
				molecule_indexes_this_domain,
				replisome_unique_indexes[replisome_parent_domain_mask]
				))

		# Get attributes of chromosomal segments on this domain
		boundary_molecule_indexes_this_domain, boundary_coordinates_this_domain, linking_numbers_this_domain = get_chromosomal_segment_attributes(
			molecule_coordinates_this_domain, molecule_indexes_this_domain,
			domain_spans_oriC, domain_spans_terC)

		# Append to existing array of attributes
		all_boundary_molecule_indexes = np.vstack([
			all_boundary_molecule_indexes, boundary_molecule_indexes_this_domain])
		all_boundary_coordinates = np.vstack((
			all_boundary_coordinates, boundary_coordinates_this_domain))
		all_segment_domain_indexes = np.concatenate((
			all_segment_domain_indexes,
			np.full(len(linking_numbers_this_domain), domain_index,
				dtype=np.int32)))
		all_linking_numbers = np.concatenate((
			all_linking_numbers, linking_numbers_this_domain))

	# Confirm total counts of all segments
	n_segments = len(all_linking_numbers)
	assert n_segments == len(active_RNAPs) + 1.5*len(replisomes) + 1

	# Add chromosomal segments
	uniqueMolCntr.objectsNew(
		'chromosomal_segment', n_segments,
		boundary_molecule_indexes=all_boundary_molecule_indexes,
		boundary_coordinates=all_boundary_coordinates,
		domain_index=all_segment_domain_indexes,
		linking_number=all_linking_numbers,
		)


def initialize_translation(bulkMolCntr, uniqueMolCntr, sim_data, randomState):
	"""
	Activate ribosomes as unique molecules, and distribute them along lengths
	of mRNAs, while decreasing counts of unactivated ribosomal subunits (30S
	and 50S).

	Ribosomes are placed randomly across the lengths of each mRNA.
	"""
	# Load translation parameters
	currentNutrients = sim_data.conditions[sim_data.condition]['nutrients']
	fracActiveRibosome = sim_data.process.translation.ribosomeFractionActiveDict[currentNutrients]
	proteinSequences = sim_data.process.translation.translation_sequences
	protein_lengths = sim_data.process.translation.monomer_data['length'].asNumber()
	translationEfficiencies = normalize(
		sim_data.process.translation.translation_efficiencies_by_monomer)
	aaWeightsIncorporated = sim_data.process.translation.translation_monomer_weights
	endWeight = sim_data.process.translation.translation_end_weight
	cistron_lengths = sim_data.process.transcription.cistron_data['length'].asNumber(units.nt)
	TU_ids = sim_data.process.transcription.rna_data['id']
	monomer_index_to_tu_indexes = sim_data.relation.monomer_index_to_tu_indexes
	monomer_index_to_cistron_index = {
		i: sim_data.process.transcription._cistron_id_to_index[monomer['cistron_id']]
		for (i, monomer) in enumerate(sim_data.process.translation.monomer_data)
		}

	# Get attributes of RNAs
	all_RNAs = uniqueMolCntr.objectsInCollection('RNA')
	TU_index_all_RNAs, length_all_RNAs, is_mRNA, is_full_transcript_all_RNAs, unique_index_all_RNAs = all_RNAs.attrs(
		'TU_index', 'transcript_length', 'is_mRNA', 'is_full_transcript', 'unique_index')
	TU_index_mRNAs = TU_index_all_RNAs[is_mRNA]
	length_mRNAs = length_all_RNAs[is_mRNA]
	is_full_transcript_mRNAs = is_full_transcript_all_RNAs[is_mRNA]
	unique_index_mRNAs = unique_index_all_RNAs[is_mRNA]

	# Calculate available template lengths of each mRNA cistron from fully
	# transcribed mRNA transcription units
	TU_index_full_mRNAs = TU_index_mRNAs[is_full_transcript_mRNAs]
	TU_counts_full_mRNAs = np.bincount(TU_index_full_mRNAs, minlength=len(TU_ids))
	cistron_counts_full_mRNAs = sim_data.process.transcription.cistron_tu_mapping_matrix.dot(
		TU_counts_full_mRNAs)
	available_cistron_lengths = np.multiply(
		cistron_counts_full_mRNAs, cistron_lengths)

	# Add available template lengths from each partially transcribed mRNAs
	TU_index_incomplete_mRNAs = TU_index_mRNAs[np.logical_not(is_full_transcript_mRNAs)]
	length_incomplete_mRNAs = length_mRNAs[np.logical_not(is_full_transcript_mRNAs)]

	TU_index_to_mRNA_lengths = {}
	for (TU_index, length) in zip(TU_index_incomplete_mRNAs, length_incomplete_mRNAs):
		TU_index_to_mRNA_lengths.setdefault(TU_index, []).append(length)

	for (TU_index, available_lengths) in TU_index_to_mRNA_lengths.items():
		cistron_indexes = sim_data.process.transcription.rna_id_to_cistron_indexes(
			TU_ids[TU_index])
		cistron_start_positions = np.array([
			sim_data.process.transcription.cistron_start_end_pos_in_tu[(cistron_index, TU_index)][0]
			for cistron_index in cistron_indexes
			])

		for length in available_lengths:
			available_cistron_lengths[cistron_indexes] += np.clip(
				length - cistron_start_positions,
				0, cistron_lengths[cistron_indexes])

	# Find number of ribosomes to activate
	ribosome30S = bulkMolCntr.countsView([sim_data.molecule_ids.s30_full_complex]).counts()[0]
	ribosome50S = bulkMolCntr.countsView([sim_data.molecule_ids.s50_full_complex]).counts()[0]
	inactiveRibosomeCount = np.minimum(ribosome30S, ribosome50S)
	n_ribosomes_to_activate = np.int64(fracActiveRibosome * inactiveRibosomeCount)

	# Add total available template lengths as weights and normalize
	protein_init_probs = normalize(
		available_cistron_lengths[sim_data.relation.cistron_to_monomer_mapping]
		* translationEfficiencies)

	# Sample a multinomial distribution of synthesis probabilities to determine
	# which types of mRNAs are initialized
	n_new_proteins = randomState.multinomial(
		n_ribosomes_to_activate, protein_init_probs)

	# Build attributes for active ribosomes
	protein_indexes = np.empty(n_ribosomes_to_activate, np.int64)
	cistron_start_positions_on_mRNA = np.empty(n_ribosomes_to_activate, np.int64)
	positions_on_mRNA_from_cistron_start_site = np.empty(n_ribosomes_to_activate, np.int64)
	mRNA_indexes = np.empty(n_ribosomes_to_activate, np.int64)
	start_index = 0
	nonzeroCount = (n_new_proteins > 0)

	for protein_index, counts in zip(
			np.arange(n_new_proteins.size)[nonzeroCount],
			n_new_proteins[nonzeroCount]):
		# Set protein index
		protein_indexes[start_index:start_index+counts] = protein_index

		# Get index of cistron corresponding to this protein
		cistron_index = monomer_index_to_cistron_index[protein_index]

		# Initialize list of available lengths for each transcript and the
		# indexes of each transcript in the list of mRNA attributes
		available_lengths = []
		attribute_indexes = []
		cistron_start_positions = []

		# Distribute ribosomes among mRNAs that produce this protein, weighted
		# by their lengths
		for TU_index in monomer_index_to_tu_indexes[protein_index]:
			attribute_indexes_this_TU = np.where(TU_index_mRNAs == TU_index)[0]
			cistron_start_position = sim_data.process.transcription.cistron_start_end_pos_in_tu[
				(cistron_index, TU_index)][0]
			available_lengths.extend(np.clip(
				length_mRNAs[attribute_indexes_this_TU] - cistron_start_position,
				0, cistron_lengths[cistron_index]))
			attribute_indexes.extend(attribute_indexes_this_TU)
			cistron_start_positions.extend(
				[cistron_start_position] * len(attribute_indexes_this_TU))

		available_lengths = np.array(available_lengths)
		attribute_indexes =	np.array(attribute_indexes)
		cistron_start_positions = np.array(cistron_start_positions)

		n_ribosomes_per_RNA = randomState.multinomial(
			counts, normalize(available_lengths))

		# Get unique indexes of each mRNA
		mRNA_indexes[start_index:start_index+counts] = np.repeat(
			unique_index_mRNAs[attribute_indexes], n_ribosomes_per_RNA)

		# Get full length of this polypeptide
		peptide_full_length = protein_lengths[protein_index]

		# Randomly place ribosomes along the length of each mRNA, capped by the
		# mRNA length expected from the full polypeptide length to prevent
		# ribosomes from overshooting full peptide lengths
		cistron_start_positions_on_mRNA[start_index:start_index+counts] = np.repeat(
			cistron_start_positions, n_ribosomes_per_RNA)
		positions_on_mRNA_from_cistron_start_site[start_index:start_index+counts] = np.floor(
			randomState.rand(counts) * np.repeat(np.minimum(available_lengths, peptide_full_length*3), n_ribosomes_per_RNA)
			)

		start_index += counts

	# Calculate the lengths of the partial polypeptide, and rescale position on
	# mRNA to be a multiple of three using this peptide length
	peptide_lengths = np.floor_divide(positions_on_mRNA_from_cistron_start_site, 3)
	positions_on_mRNA = cistron_start_positions_on_mRNA + 3*peptide_lengths

	# Update masses of partially translated proteins
	sequences = proteinSequences[protein_indexes]
	mass_increase_protein = computeMassIncrease(
		sequences, peptide_lengths, aaWeightsIncorporated)

	# Add end weight
	mass_increase_protein[peptide_lengths != 0] += endWeight

	# Add active ribosomes
	uniqueMolCntr.objectsNew(
		'active_ribosome', n_ribosomes_to_activate,
		protein_index=protein_indexes,
		peptide_length=peptide_lengths,
		mRNA_index=mRNA_indexes,
		pos_on_mRNA=positions_on_mRNA,
		massDiff_protein=mass_increase_protein,
		)

	# Decrease counts of free 30S and 50S ribosomal subunits
	bulkMolCntr.countsIs(ribosome30S - n_ribosomes_to_activate,
		[sim_data.molecule_ids.s30_full_complex])
	bulkMolCntr.countsIs(ribosome50S - n_ribosomes_to_activate,
		[sim_data.molecule_ids.s50_full_complex])


def determine_chromosome_state(C, D, tau, replichore_length, n_max_replisomes,
		place_holder):
	"""
	Calculates the attributes of oriC's, replisomes, and chromosome domains on
	the chromosomes at the beginning of the cell cycle.

	Inputs
	--------
	- C: the C period of the cell, the length of time between replication
	initiation and replication termination.
	- D: the D period of the cell, the length of time between completing
	replication of the chromosome and division of the cell.
	- tau: the doubling time of the cell
	- n_max_replisomes: the maximum number of replisomes that can be formed
	given the initial counts of replisome subunits
	- replichore_length: the amount of DNA to be replicated per fork, usually
	half of the genome, in base-pairs
	- place_holder: placeholder value for chromosome domains without child
	domains

	Returns
	--------
	oric_state: dictionary of attributes for the oriC molecules with the
	following keys.
	- domain_index: a vector of integers indicating which chromosome domain the
	oriC sequence belongs to.

	replisome_state: dictionary of attributes for the replisome molecules
	with the following keys.
	- coordinates: a vector of integers that indicates where the replisomes
	are located on the chromosome relative to the origin, in base pairs.
	- right_replichore: a vector of boolean values that indicates whether the
	replisome is on the right replichore (True) or the left replichore (False).
	- domain_index: a vector of integers indicating which chromosome domain the
	replisomes belong to. The index of the "mother" domain of the replication
	fork is assigned to the replisome.

	domain_state: dictionary of attributes for the chromosome domains with the
	following keys.
	- domain_index: the indexes of the domains.
	- child_domains: the (n_domain X 2) array of the domain indexes of the two
	children domains that are connected on the oriC side with the given domain.
	"""

	# All inputs must be positive numbers
	assert C.asNumber(units.min) >= 0, "C value can't be negative."
	assert D.asNumber(units.min) >= 0, "D value can't be negative."
	assert tau.asNumber(units.min) >= 0, "tau value can't be negative."
	assert replichore_length.asNumber(units.nt) >= 0, "replichore_length value can't be negative."

	# Require that D is shorter than tau - time between completing DNA
	# replication and cell division must be shorter than the time between two
	# cell divisions.
	assert D.asNumber(units.min) < tau.asNumber(units.min), "The D period must be shorter than the doubling time tau."

	# Calculate the maximum number of replication rounds given the maximum
	# count of replisomes
	n_max_rounds = int(np.log2(n_max_replisomes/2 + 1))

	# Calculate the number of active replication rounds
	n_rounds = min(n_max_rounds,
		int(np.floor(
		(C.asNumber(units.min) + D.asNumber(units.min))/tau.asNumber(units.min)
		)))

	# Initialize arrays for replisomes
	n_replisomes = 2*(2**n_rounds - 1)
	coordinates = np.zeros(n_replisomes, dtype=np.int64)
	right_replichore_replisome = np.zeros(n_replisomes, dtype=bool)
	domain_index_replisome = np.zeros(n_replisomes, dtype=np.int32)

	# Initialize child domain array for chromosome domains
	n_domains = 2**(n_rounds + 1) - 1
	child_domains = np.full((n_domains, 2), place_holder, dtype=np.int32)

	# Set domain_index attribute of oriC's and chromosome domains
	domain_index_oric = np.arange(2**n_rounds - 1, 2**(n_rounds + 1) - 1, dtype=np.int32)
	domain_index_domains = np.arange(0, n_domains, dtype=np.int32)

	def n_events_before_this_round(round_idx):
		"""
		Calculates the number of replication events that happen before the
		replication round index given as an argument. Since 2**i events happen
		at each round i = 0, 1, ..., the sum of the number of events before
		round j is 2**j - 1.
		"""
		return 2**round_idx - 1

	# Loop through active replication rounds, starting from the oldest round.
	# If n_round = 0 skip loop entirely - no active replication round.
	for round_idx in np.arange(n_rounds):
		# Determine at which location (base) of the chromosome the replication
		# forks should be initialized to
		rel_location = 1.0 - (((round_idx + 1.0)*tau - D)/C)
		rel_location = units.strip_empty_units(rel_location)
		fork_location = np.floor(rel_location*(
			replichore_length.asNumber(units.nt)))

		# Add 2^n initiation events per round. A single initiation event
		# generates two replication forks.
		n_events_this_round = 2**round_idx

		# Set attributes of replisomes for this replication round
		coordinates[
			2*n_events_before_this_round(round_idx):
			2*n_events_before_this_round(round_idx + 1)
			] = np.tile(np.array([fork_location, -fork_location]), n_events_this_round)

		right_replichore_replisome[
			2*n_events_before_this_round(round_idx):
			2*n_events_before_this_round(round_idx + 1)
			] = np.tile(np.array([True, False]), n_events_this_round)

		for i, domain_index in enumerate(
				np.arange(n_events_before_this_round(round_idx),
					n_events_before_this_round(round_idx+1))):
			domain_index_replisome[
				2*n_events_before_this_round(round_idx) + 2*i:
				2*n_events_before_this_round(round_idx) + 2*(i+1)
				] = np.repeat(domain_index, 2)

		# Set attributes of chromosome domains for this replication round
		for i, domain_index in enumerate(
				np.arange(n_events_before_this_round(round_idx + 1),
					n_events_before_this_round(round_idx + 2), 2)):
			child_domains[
				n_events_before_this_round(round_idx) + i, :
				] = np.array([domain_index, domain_index + 1])

	# Convert to numpy arrays and wrap into dictionaries
	oric_state = {"domain_index": domain_index_oric}

	replisome_state = {
		"coordinates": coordinates,
		"right_replichore": right_replichore_replisome,
		"domain_index": domain_index_replisome}

	domain_state = {
		"child_domains": child_domains,
		"domain_index": domain_index_domains}

	return oric_state, replisome_state, domain_state


def rescale_initiation_probs(init_probs, TU_index, fixed_synth_probs,
		fixed_TU_indexes):
	"""
	Rescales the initiation probabilities of each promoter such that the total
	synthesis probabilities of certain types of RNAs are fixed to a
	predetermined value. For instance, if there are two copies of promoters for
	RNA A, whose synthesis probability should be fixed to 0.1, each promoter is
	given an initiation probability of 0.05.
	"""
	for rna_idx, synth_prob in zip(fixed_TU_indexes, fixed_synth_probs):
		fixed_rna_mask = (TU_index == rna_idx)
		init_probs[fixed_rna_mask] = synth_prob / fixed_rna_mask.sum()
