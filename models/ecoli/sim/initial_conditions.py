"""
TODO: document math
TODO: raise/warn if physiological metabolite concentration targets appear to be smaller than what
  is needed at this time step size
"""

from __future__ import absolute_import, division, print_function

from itertools import izip
import scipy.sparse

import numpy as np

from wholecell.utils.fitting import normalize, countsFromMassAndExpression, masses_and_counts_for_homeostatic_target
from wholecell.utils.polymerize import computeMassIncrease
from wholecell.utils import units
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices
from wholecell.utils.random import stochasticRound

from models.ecoli.processes.polypeptide_elongation import SteadyStateElongationModel
from wholecell.containers.unique_objects_container import Access
from wholecell.sim.divide_cell import load_inherited_state


RAND_MAX = 2**31

def calcInitialConditions(sim, sim_data):
	'''Calculate the initial conditions for a new cell without inherited state
	from a parent cell.

	Params:
		sim: The simulation object, which must not have an _inheritedStatePath,
			otherwise the simulation should use setDaughterInitialConditions()
		sim_data: to initialize state from
	'''
	assert sim._inheritedStatePath is None
	randomState = sim.randomState

	massCoeff = 1.0
	if sim._massDistribution:
		massCoeff = randomState.normal(loc = 1.0, scale = 0.1)

	bulkMolCntr = sim.internal_states['BulkMolecules'].container
	uniqueMolCntr = sim.internal_states["UniqueMolecules"].container

	# Set up states
	initializeBulkMolecules(bulkMolCntr, sim_data, sim.external_states['Environment'].current_media_id,
		randomState, massCoeff, sim._ppgpp_regulation)
	initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

	# Must be called after unique and bulk molecules are initialized to get
	# concentrations for ribosomes, tRNA, synthetases etc from cell volume
	if sim._trna_charging:
		elongation_model = SteadyStateElongationModel(sim_data, sim.processes['PolypeptideElongation'])
		initialize_trna_charging(sim_data, sim.internal_states, elongation_model.calculate_trna_charging)

def initializeBulkMolecules(bulkMolCntr, sim_data, current_media_id, randomState, massCoeff, ppgpp_regulation):

	# Set protein counts from expression
	initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff, ppgpp_regulation)

	# Set RNA counts from expression
	initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff, ppgpp_regulation)

	# Set other biomass components
	initializeSmallMolecules(bulkMolCntr, sim_data, current_media_id, randomState, massCoeff)

	# Form complexes
	initializeComplexation(bulkMolCntr, sim_data, randomState)

def initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, randomState):
	# Initialize counts of full chromosomes
	initializeFullChromosome(bulkMolCntr, uniqueMolCntr, sim_data)

	# Initialize unique molecules relevant to replication
	initializeReplication(bulkMolCntr, uniqueMolCntr, sim_data)

	# Initialize bound transcription factors
	initialize_transcription_factors(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

	# Initialize active RNAPs and unique molecule representations of RNAs
	initialize_transcription(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

	# Initialize activate ribosomes
	initialize_translation(bulkMolCntr, uniqueMolCntr, sim_data, randomState)

def initialize_trna_charging(sim_data, states, calc_charging):
	'''
	Initializes charged tRNA from uncharged tRNA and amino acids

	Inputs:
		sim_data (SimulationDataEcoli object)
		states (dict with internal_state objects as values) - internal states of sim
		calc_charging (function) - function to calculate charging of tRNA

	Notes:
		Does not adjust for mass of amino acids on charged tRNA (~0.01% of cell mass)
	'''

	# Calculate cell volume for concentrations
	mass = 0
	for state in states.values():
		state.calculateMass()
		mass += np.sum(state.mass())
	cell_volume = units.fg * mass / sim_data.constants.cellDensity
	counts_to_molar = 1 / (sim_data.constants.nAvogadro * cell_volume)

	# Get molecule views and concentrations
	transcription = sim_data.process.transcription
	aa_from_synthetase = transcription.aa_from_synthetase
	aa_from_trna = transcription.aa_from_trna
	bulk_molecules = states['BulkMolecules'].container
	synthetases = bulk_molecules.countsView(transcription.synthetase_names)
	uncharged_trna = bulk_molecules.countsView(transcription.rnaData['id'][transcription.rnaData['isTRna']])
	charged_trna = bulk_molecules.countsView(transcription.charged_trna_names)
	aas = bulk_molecules.countsView(sim_data.moleculeGroups.aaIDs)
	ribosome_counts = states['UniqueMolecules'].container.counts(['active_ribosome'])

	synthetase_conc = counts_to_molar * np.dot(aa_from_synthetase, synthetases.counts())
	uncharged_trna_conc = counts_to_molar * np.dot(aa_from_trna, uncharged_trna.counts())
	charged_trna_conc = counts_to_molar * np.dot(aa_from_trna, charged_trna.counts())
	aa_conc = counts_to_molar * aas.counts()
	ribosome_conc = counts_to_molar * ribosome_counts

	# Estimate fraction of amino acids from sequences, excluding first index for padding of -1
	_, aas_in_sequences = np.unique(sim_data.process.translation.translationSequences, return_counts=True)
	f = aas_in_sequences[1:] / np.sum(aas_in_sequences[1:])

	# Estimate initial charging state
	fraction_charged, _ = calc_charging(synthetase_conc, uncharged_trna_conc, charged_trna_conc, aa_conc, ribosome_conc, f)

	# Update counts of tRNA to match charging
	total_trna_counts = uncharged_trna.counts() + charged_trna.counts()
	charged_trna_counts = np.round(total_trna_counts * np.dot(fraction_charged, aa_from_trna))
	uncharged_trna_counts = total_trna_counts - charged_trna_counts
	charged_trna.countsIs(charged_trna_counts)
	uncharged_trna.countsIs(uncharged_trna_counts)

def initializeProteinMonomers(bulkMolCntr, sim_data, randomState, massCoeff, ppgpp_regulation):

	monomersView = bulkMolCntr.countsView(sim_data.process.translation.monomerData["id"])
	monomerMass = massCoeff * sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])["proteinMass"] / sim_data.mass.avgCellToInitialCellConvFactor
	# TODO: unify this logic with the parca so it doesn't fall out of step
	# again (look at the calcProteinCounts function)

	if ppgpp_regulation:
		ppgpp = sim_data.growthRateParameters.getppGppConc(sim_data.conditionToDoublingTime[sim_data.condition])
		rnaExpression = sim_data.process.transcription.expression_from_ppgpp(ppgpp)
	else:
		rnaExpression = sim_data.process.transcription.rnaExpression[sim_data.condition]

	monomerExpression = normalize(
		rnaExpression[sim_data.relation.rnaIndexToMonomerMapping] *
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


def initializeRNA(bulkMolCntr, sim_data, randomState, massCoeff, ppgpp_regulation):
	"""
	Initializes counts of RNAs in the bulk molecule container using RNA
	expression data. mRNA counts are also initialized here, but is later reset
	to zero when the representations for mRNAs are moved to the unique molecule
	container.
	"""
	rnaView = bulkMolCntr.countsView(sim_data.process.transcription.rnaData["id"])
	rnaMass = massCoeff * sim_data.mass.getFractionMass(sim_data.conditionToDoublingTime[sim_data.condition])["rnaMass"] / sim_data.mass.avgCellToInitialCellConvFactor

	if ppgpp_regulation:
		ppgpp = sim_data.growthRateParameters.getppGppConc(sim_data.conditionToDoublingTime[sim_data.condition])
		rnaExpression = sim_data.process.transcription.expression_from_ppgpp(ppgpp)
	else:
		rnaExpression = normalize(sim_data.process.transcription.rnaExpression[sim_data.condition])

	nRnas = countsFromMassAndExpression(
		rnaMass.asNumber(units.g),
		sim_data.process.transcription.rnaData["mw"].asNumber(units.g/units.mol),
		rnaExpression,
		sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		)

	# ID Groups of rRNAs
	idx_16Srrna = np.where(sim_data.process.transcription.rnaData['isRRna16S'])[0]
	idx_23Srrna = np.where(sim_data.process.transcription.rnaData['isRRna23S'])[0]
	idx_5Srrna = np.where(sim_data.process.transcription.rnaData['isRRna5S'])[0]

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
def initializeSmallMolecules(bulkMolCntr, sim_data, current_media_id, randomState, massCoeff):
	doubling_time = sim_data.conditionToDoublingTime[sim_data.condition]
	avgCellFractionMass = sim_data.mass.getFractionMass(doubling_time)

	mass = massCoeff * (avgCellFractionMass["proteinMass"] + avgCellFractionMass["rnaMass"] + avgCellFractionMass["dnaMass"]) / sim_data.mass.avgCellToInitialCellConvFactor

	concDict = sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
		current_media_id
		)
	concDict.update(sim_data.mass.getBiomassAsConcentrations(doubling_time))
	concDict[sim_data.moleculeIds.ppGpp] = sim_data.growthRateParameters.getppGppConc(doubling_time)
	moleculeIds = sorted(concDict)
	moleculeConcentrations = (units.mol / units.L) * np.array([concDict[key].asNumber(units.mol / units.L) for key in moleculeIds])

	massesToAdd, countsToAdd = masses_and_counts_for_homeostatic_target(
		mass,
		moleculeConcentrations,
		sim_data.getter.getMass(moleculeIds),
		sim_data.constants.cellDensity,
		sim_data.constants.nAvogadro
		)

	bulkMolCntr.countsIs(
		countsToAdd,
		moleculeIds
		)

def initializeComplexation(bulkMolCntr, sim_data, randomState):
	moleculeNames = sim_data.process.complexation.moleculeNames
	moleculeView = bulkMolCntr.countsView(moleculeNames)

	# save rnase counts for after complexation, need to be monomers
	rnases = sim_data.process.rna_decay.endoRnaseIds
	rnaseCounts = bulkMolCntr.countsView(rnases).counts()
	bulkMolCntr.countsIs(0, rnases)

	stoichMatrix = sim_data.process.complexation.stoichMatrix().astype(np.int64, order='F')
	prebuiltMatrices = mccBuildMatrices(stoichMatrix)

	# form complexes until no new complexes form (some complexes are complexes of complexes)
	while True:
		moleculeCounts = moleculeView.counts()
		updatedMoleculeCounts, complexationEvents = mccFormComplexesWithPrebuiltMatrices(
			moleculeCounts,
			randomState.randint(1000),
			stoichMatrix,
			*prebuiltMatrices)

		bulkMolCntr.countsIs(
			updatedMoleculeCounts,
			moleculeNames)

		if np.any(updatedMoleculeCounts < 0):
			raise ValueError('Negative counts after complexation')

		if not np.any(moleculeCounts - updatedMoleculeCounts):
			break

	bulkMolCntr.countsIs(rnaseCounts, rnases)


def initializeFullChromosome(bulkMolCntr, uniqueMolCntr, sim_data):
	"""
	Initializes the counts of full chromosomes to one. The division_time of
	this initial chromosome is set to be zero for consistency.
	"""
	uniqueMolCntr.objectsNew(
		'full_chromosome', 1,
		division_time=0.0,
		has_induced_division=True,
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
	C = sim_data.growthRateParameters.c_period.asUnit(units.min)
	D = sim_data.growthRateParameters.d_period.asUnit(units.min)
	tau = sim_data.conditionToDoublingTime[sim_data.condition].asUnit(units.min)

	# Calculate length of replichore
	genome_length = sim_data.process.replication.genome_length
	replichore_length = np.ceil(0.5*genome_length) * units.nt

	# Calculate the maximum number of replisomes that could be formed with
	# the existing counts of replisome subunits
	n_max_replisomes = np.min(np.concatenate(
		(bulkMolCntr.counts(sim_data.moleculeGroups.replisome_trimer_subunits)//3,
		bulkMolCntr.counts(sim_data.moleculeGroups.replisome_monomer_subunits))))

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
			sim_data.process.replication.replicationMonomerWeights.asNumber(units.fg))

		# Add active replisomes as unique molecules and set attributes
		uniqueMolCntr.objectsNew(
			'active_replisome', n_replisome,
			coordinates=replisome_state["coordinates"],
			right_replichore=replisome_state["right_replichore"],
			domain_index=replisome_state["domain_index"],
			massDiff_DNA=mass_increase_dna[0::2] + mass_increase_dna[1::2])

		# Remove replisome subunits from bulk molecules
		bulkMolCntr.countsDec(3*n_replisome, sim_data.moleculeGroups.replisome_trimer_subunits)
		bulkMolCntr.countsDec(n_replisome, sim_data.moleculeGroups.replisome_monomer_subunits)

	# Get coordinates of all promoters and DnaA boxes
	all_promoter_coordinates = sim_data.process.transcription.rnaData[
		"replicationCoordinate"]
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
					motif_mask = np.ones_like(all_motif_coordinates, dtype=np.bool)

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
		bound_TF=np.zeros((n_promoter, n_tf), dtype=np.bool))

	# Add DnaA boxes as unique molecules and set attributes
	n_DnaA_box = len(DnaA_box_coordinates)

	uniqueMolCntr.objectsNew(
		'DnaA_box', n_DnaA_box,
		coordinates=np.array(DnaA_box_coordinates),
		domain_index=np.array(DnaA_box_domain_index),
		DnaA_bound=np.zeros(n_DnaA_box, dtype=np.bool)
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
	tfToTfType = sim_data.process.transcription_regulation.tfToTfType
	pPromoterBoundTF = sim_data.process.transcription_regulation.pPromoterBoundTF

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
			if tf == sim_data.process.transcription_regulation.activeToBound[tf]:
				inactive_tf_view[tf] = bulkMolCntr.countsView([
					sim_data.process.equilibrium.getUnbound(tf + "[c]")
					])
			else:
				inactive_tf_view[tf] = bulkMolCntr.countsView([
					sim_data.process.transcription_regulation.activeToBound[tf] + "[c]"
					])
		elif tfToTfType[tf] == "2CS":
			inactive_tf_view[tf] = bulkMolCntr.countsView([
				sim_data.process.two_component_system.activeToInactiveTF[tf + "[c]"]
				])

	# Get masses of active transcription factors
	bulk_molecule_ids = sim_data.internal_state.bulkMolecules.bulkData["id"]
	tf_indexes = [np.where(bulk_molecule_ids == tf_id + "[c]")[0][0]
		for tf_id in tf_ids]
	active_tf_masses = (sim_data.internal_state.bulkMolecules.bulkData["mass"][
		tf_indexes] / sim_data.constants.nAvogadro).asNumber(units.fg)

	# Get attributes of promoters
	promoters = uniqueMolCntr.objectsInCollection(
		'promoter', access=[Access.EDIT])
	TU_index = promoters.attr("TU_index")

	# Initialize bound_TF array
	bound_TF = np.zeros((len(promoters), len(tf_ids)), dtype=np.bool)

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
			randomState, n_available_promoters * pPromoterBound))

		bound_locs = np.zeros(n_available_promoters, dtype=np.bool)
		if n_to_bind > 0:
			# Determine randomly which DNA targets to bind based on which of
			# the following is more limiting:
			# number of promoter sites to bind, or number of active
			# transcription factors
			bound_locs[
				randomState.choice(
					n_available_promoters,
					size=np.min((n_to_bind, active_tf_view[tf_id].counts())),
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


def initialize_transcription(bulkMolCntr, uniqueMolCntr, sim_data, randomState):
	"""
	Activate RNA polymerases as unique molecules, and distribute them along
	length of genes, while decreasing counts of unactivated RNA polymerases
	(APORNAP-CPLX[c]). Also initialize unique molecule representations of fully
	transcribed mRNAs and partially transcribed RNAs, using counts of mRNAs
	initialized as bulk molecules, and the attributes of initialized RNA
	polymerases. The counts of full mRNAs represented as bulk molecules are
	reset to zero.

	RNA polymerases are placed randomly across the length of each gene, with
	the synthesis probabilities for each gene determining the number of RNA
	polymerases placed at each gene.
	"""
	# Load parameters
	rnaLengths = sim_data.process.transcription.rnaData['length'].asNumber()
	rna_masses = (sim_data.process.transcription.rnaData['mw']/sim_data.constants.nAvogadro).asNumber(units.fg)
	current_media_id = sim_data.conditions[sim_data.condition]['nutrients']
	fracActiveRnap = sim_data.process.transcription.rnapFractionActiveDict[current_media_id]
	inactive_RNAP_counts = bulkMolCntr.countsView([sim_data.moleculeIds.rnapFull]).counts()[0]
	rnaSequences = sim_data.process.transcription.transcriptionSequences
	ntWeights = sim_data.process.transcription.transcriptionMonomerWeights
	endWeight = sim_data.process.transcription.transcriptionEndWeight
	replichore_lengths = sim_data.process.replication.replichore_lengths
	chromosome_length = replichore_lengths.sum()

	# Number of rnaPoly to activate
	n_RNAPs_to_activate = np.int64(fracActiveRnap * inactive_RNAP_counts)

	# Parameters for rnaSynthProb
	basal_prob = sim_data.process.transcription_regulation.basal_prob
	n_TUs = len(basal_prob)
	delta_prob = sim_data.process.transcription_regulation.delta_prob
	delta_prob_matrix = scipy.sparse.csr_matrix(
		(delta_prob['deltaV'], (delta_prob['deltaI'], delta_prob['deltaJ'])),
		shape=delta_prob['shape']).toarray()

	# Get attributes of promoters
	promoters = uniqueMolCntr.objectsInCollection("promoter")
	n_promoters = len(promoters)
	TU_index, bound_TF, domain_index_promoters = promoters.attrs(
		"TU_index", "bound_TF", "domain_index")

	# Construct matrix that maps promoters to transcription units
	TU_to_promoter = scipy.sparse.csr_matrix(
		(np.ones(n_promoters), (TU_index, np.arange(n_promoters))),
		shape=(n_TUs, n_promoters))

	# Synthesis probabilities for different categories of genes
	rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction
	rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
	rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase

	# Get coordinates and transcription directions of transcription units
	replication_coordinate = sim_data.process.transcription.rnaData[
		"replicationCoordinate"]
	transcription_direction = sim_data.process.transcription.rnaData[
		"direction"]

	# Determine changes from genetic perturbations
	genetic_perturbations = {}
	perturbations = getattr(sim_data, 'genetic_perturbations', {})

	if len(perturbations) > 0:
		probability_indexes = [
			(index, sim_data.genetic_perturbations[rna_data['id']])
				for index, rna_data in enumerate(sim_data.process.transcription.rnaData)
				if rna_data['id'] in sim_data.genetic_perturbations]

		genetic_perturbations = {
			'fixedRnaIdxs': map(lambda pair: pair[0], probability_indexes),
			'fixedSynthProbs': map(lambda pair: pair[1], probability_indexes)}

	# If initiationShuffleIdxs does not exist, set value to None
	shuffleIdxs = getattr(sim_data.process.transcription, 'initiationShuffleIdxs', None)

	# ID Groups
	idx_rRNA = np.where(sim_data.process.transcription.rnaData['isRRna'])[0]
	idx_mRNA = np.where(sim_data.process.transcription.rnaData["isMRna"])[0]
	idx_tRNA = np.where(sim_data.process.transcription.rnaData["isTRna"])[0]
	idx_rprotein = np.where(sim_data.process.transcription.rnaData['isRProtein'])[0]
	idx_rnap = np.where(sim_data.process.transcription.rnaData['isRnap'])[0]

	# Calculate probabilities of the RNAP binding to the promoters
	promoter_init_probs = (basal_prob[TU_index] +
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

	# Randomly advance RNAPs along the gene
	# TODO (Eran): make sure there aren't any RNAPs at same location on same gene
	updated_lengths = np.array(
		randomState.rand(n_RNAPs_to_activate) * rnaLengths[TU_index_partial_RNAs],
		dtype=np.int)

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
		[sim_data.moleculeIds.rnapFull])

	# Add partially transcribed RNAs
	uniqueMolCntr.objectsNew(
		'RNA', n_RNAPs_to_activate,
		TU_index=TU_index_partial_RNAs,
		transcript_length=updated_lengths,
		is_mRNA=is_mRNA_partial_RNAs,
		is_full_transcript=np.zeros(n_RNAPs_to_activate, dtype=np.bool),
		can_translate=is_mRNA_partial_RNAs,
		RNAP_index=RNAP_indexes,
		massDiff_RNA=added_RNA_mass,
		massDiff_mRNA=added_mRNA_mass)

	# Get counts of mRNAs initialized as bulk molecules
	mRNA_ids = sim_data.process.transcription.rnaData["id"][
		sim_data.process.transcription.rnaData["isMRna"]]
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
		is_mRNA=np.ones_like(TU_index_full_mRNAs, dtype=np.bool),
		is_full_transcript=np.ones_like(TU_index_full_mRNAs, dtype=np.bool),
		can_translate=np.ones_like(TU_index_full_mRNAs, dtype=np.bool),
		RNAP_index=np.full(TU_index_full_mRNAs.shape, -1, dtype=np.int64),
		massDiff_mRNA=rna_masses[TU_index_full_mRNAs])

	# Reset counts of bulk mRNAs to zero
	mRNA_view.countsIs(0)


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
	proteinSequences = sim_data.process.translation.translationSequences
	translationEfficiencies = normalize(
		sim_data.process.translation.translationEfficienciesByMonomer)
	aaWeightsIncorporated = sim_data.process.translation.translationMonomerWeights
	endWeight = sim_data.process.translation.translationEndWeight

	# Get attributes of RNAs
	all_RNAs = uniqueMolCntr.objectsInCollection('RNA')
	TU_index_all_RNAs, length_all_RNAs, is_mRNA, unique_index_all_RNAs = all_RNAs.attrs(
		'TU_index', 'transcript_length', 'is_mRNA', 'unique_index')
	TU_index_mRNAs = TU_index_all_RNAs[is_mRNA]
	length_mRNAs = length_all_RNAs[is_mRNA]
	unique_index_mRNAs = unique_index_all_RNAs[is_mRNA]

	# Get conversion matrix between transcription units and mRNA templates
	# for each monomer
	# TODO (ggsun): This should be modified when transcription unit structures
	# 	are incorporated - the matrix would not be one-on-one, and the lengths
	# 	of each transcript will affect the availabilities of mRNAs on that
	# 	transcript.
	all_TU_ids = sim_data.process.transcription.rnaData['id']
	all_mRNA_ids = sim_data.process.translation.monomerData['rnaId']
	TU_counts_to_mRNA_counts = np.zeros(
		(len(all_mRNA_ids), len(all_TU_ids)), dtype=np.int64)

	TU_id_to_index = {TU_id: i for i, TU_id in enumerate(all_TU_ids)}
	protein_index_to_TU_index = {}
	for i, mRNA_id in enumerate(all_mRNA_ids):
		TU_counts_to_mRNA_counts[i, TU_id_to_index[mRNA_id]] = 1
		protein_index_to_TU_index[i] = TU_id_to_index[mRNA_id]

	# Calculate available template lengths of each mRNA
	TU_total_length = np.zeros(len(all_TU_ids), dtype=np.int64)
	for index, length in zip(TU_index_mRNAs, length_mRNAs):
		TU_total_length[index] += length

	mRNA_total_length = TU_counts_to_mRNA_counts.dot(TU_total_length)

	# Find number of ribosomes to activate
	ribosome30S = bulkMolCntr.countsView([sim_data.moleculeIds.s30_fullComplex]).counts()[0]
	ribosome50S = bulkMolCntr.countsView([sim_data.moleculeIds.s50_fullComplex]).counts()[0]
	inactiveRibosomeCount = np.minimum(ribosome30S, ribosome50S)
	n_ribosomes_to_activate = np.int64(fracActiveRibosome * inactiveRibosomeCount)

	# Add total available template lengths as weights and normalize
	protein_init_probs = normalize(mRNA_total_length*translationEfficiencies)

	# Sample a multinomial distribution of synthesis probabilities to determine
	# which types of mRNAs are initialized
	n_new_proteins = randomState.multinomial(
		n_ribosomes_to_activate, protein_init_probs)

	# Build attributes for active ribosomes
	protein_indexes = np.empty(n_ribosomes_to_activate, np.int64)
	positions_on_mRNA = np.empty(n_ribosomes_to_activate, np.int64)
	mRNA_indexes = np.empty(n_ribosomes_to_activate, np.int64)
	start_index = 0
	nonzeroCount = (n_new_proteins > 0)

	for protein_index, counts in zip(
			np.arange(n_new_proteins.size)[nonzeroCount],
			n_new_proteins[nonzeroCount]):
		# Set protein index
		protein_indexes[start_index:start_index+counts] = protein_index

		# Distribute ribosomes among mRNAs that produce this protein, weighted
		# by their lengths
		mask = (TU_index_mRNAs == protein_index_to_TU_index[protein_index])
		lengths = length_mRNAs[mask]
		n_ribosomes_per_RNA = randomState.multinomial(
			counts, normalize(lengths))

		# Get unique indexes of each mRNA
		mRNA_indexes[start_index:start_index+counts] = np.repeat(
			unique_index_mRNAs[mask], n_ribosomes_per_RNA)

		# Randomly place ribosomes along the length of each mRNA
		positions_on_mRNA[start_index:start_index+counts] = np.floor(
			randomState.rand(counts)*np.repeat(lengths, n_ribosomes_per_RNA))

		start_index += counts

	# Calculate the lengths of the partial polypeptide, and rescale position on
	# mRNA to be a multiple of three using this peptide length
	peptide_lengths = np.floor_divide(positions_on_mRNA, 3)
	positions_on_mRNA = 3*peptide_lengths

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
		[sim_data.moleculeIds.s30_fullComplex])
	bulkMolCntr.countsIs(ribosome50S - n_ribosomes_to_activate,
		[sim_data.moleculeIds.s50_fullComplex])


def setDaughterInitialConditions(sim, sim_data):
	'''Calculate the initial conditions for a new cell from state inherited
	from a parent cell, stored in a file.

	Params:
		sim: The simulation object, which must have an _inheritedStatePath to
			name the file directory to load the state from, otherwise the
			simulation should use calcInitialConditions()
		sim_data: Unused argument needed to conform to an informal interface
	'''
	inherited_state_path = sim._inheritedStatePath
	assert inherited_state_path is not None

	inherited_state = load_inherited_state(inherited_state_path)

	sim._isDead = inherited_state['is_dead']

	elng_rate_factor = inherited_state['elng_rate_factor']
	if sim._growthRateNoise:
		sim.processes["PolypeptideElongation"].elngRateFactor = elng_rate_factor

	sim.internal_states["BulkMolecules"].loadSnapshot(inherited_state['bulk_molecules'])
	sim.internal_states["UniqueMolecules"].loadSnapshot(inherited_state['unique_molecules'])

	sim._initialTime = inherited_state['initial_time']


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
	right_replichore_replisome = np.zeros(n_replisomes, dtype=np.bool)
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
		rel_location = units.convertNoUnitToNumber(rel_location)
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
	for rna_idx, synth_prob in izip(fixed_TU_indexes, fixed_synth_probs):
		fixed_rna_mask = (TU_index == rna_idx)
		init_probs[fixed_rna_mask] = synth_prob / fixed_rna_mask.sum()
