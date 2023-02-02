import numpy as np

from models.ecoli.processes.polypeptide_elongation import calculate_steady_state_trna_charging, get_charging_params
import reconstruction.ecoli.initialization as init
from wholecell.sim.divide_cell import load_inherited_state


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
	media_id = sim.external_states['Environment'].current_media_id
	import_molecules = sim.external_states['Environment'].get_import_molecules()

	# Set up states
	init.initializeBulkMolecules(bulkMolCntr, sim_data, media_id, import_molecules,
		randomState, massCoeff, sim._ppgpp_regulation, sim._trna_attenuation)
	cell_mass = init.calculate_cell_mass(sim.internal_states)
	init.initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data, cell_mass,
		randomState, sim._superhelical_density, sim._ppgpp_regulation, sim._trna_attenuation, sim._kinetic_trna_charging)

	# Must be called after unique and bulk molecules are initialized to get
	# concentrations for ribosomes, tRNA, synthetases etc from cell volume
	if sim._kinetic_trna_charging:
		initialize_kinetic_trna_charging(sim_data.condition, sim.internal_states['BulkMolecules'].container, sim_data, sim.internal_states)
	elif sim._steady_state_trna_charging:
		initialize_steady_state_trna_charging(sim_data, sim.internal_states, sim._variable_elongation_translation)

	# Adjust small molecule concentrations again after other mass adjustments
	# for more stable metabolism solution at beginning of sims
	init.set_small_molecule_counts(bulkMolCntr, sim_data, media_id, import_molecules,
		massCoeff, cell_mass=init.calculate_cell_mass(sim.internal_states))

def initialize_steady_state_trna_charging(sim_data, states, variable_elongation):
	'''
	Initializes charged tRNA from uncharged tRNA and amino acids

	Inputs:
		sim_data (SimulationDataEcoli object)
		states (dict with internal_state objects as values) - internal states of sim
		variable_elongation (bool) - if True, the max elongation rate is set to be
			higher in the simulation

	Notes:
		Does not adjust for mass of amino acids on charged tRNA (~0.01% of cell mass)
	'''

	# Calculate cell volume for concentrations
	cell_volume = init.calculate_cell_mass(states) / sim_data.constants.cell_density
	counts_to_molar = 1 / (sim_data.constants.n_avogadro * cell_volume)

	# Get molecule views and concentrations
	transcription = sim_data.process.transcription
	aa_from_synthetase = transcription.aa_from_synthetase
	aa_from_trna = transcription.aa_from_trna
	bulk_molecules = states['BulkMolecules'].container
	synthetases = bulk_molecules.countsView(transcription.synthetase_names)
	uncharged_trna = bulk_molecules.countsView(transcription.rna_data['id'][transcription.rna_data['is_tRNA']])
	charged_trna = bulk_molecules.countsView(transcription.charged_trna_names)
	aas = bulk_molecules.countsView(sim_data.molecule_groups.amino_acids)
	ribosome_counts = states['UniqueMolecules'].container.counts(['active_ribosome'])

	synthetase_conc = counts_to_molar * np.dot(aa_from_synthetase, synthetases.counts())
	uncharged_trna_conc = counts_to_molar * np.dot(aa_from_trna, uncharged_trna.counts())
	charged_trna_conc = counts_to_molar * np.dot(aa_from_trna, charged_trna.counts())
	aa_conc = counts_to_molar * aas.counts()
	ribosome_conc = counts_to_molar * ribosome_counts

	# Estimate fraction of amino acids from sequences, excluding first index for padding of -1
	_, aas_in_sequences = np.unique(sim_data.process.translation.translation_sequences, return_counts=True)
	f = aas_in_sequences[1:] / np.sum(aas_in_sequences[1:])

	# Estimate initial charging state
	charging_params = get_charging_params(sim_data, variable_elongation=variable_elongation)
	fraction_charged, *_ = calculate_steady_state_trna_charging(synthetase_conc, uncharged_trna_conc,
		charged_trna_conc, aa_conc, ribosome_conc, f, charging_params)

	# Update counts of tRNA to match charging
	total_trna_counts = uncharged_trna.counts() + charged_trna.counts()
	charged_trna_counts = np.round(total_trna_counts * np.dot(fraction_charged, aa_from_trna))
	uncharged_trna_counts = total_trna_counts - charged_trna_counts
	charged_trna.countsIs(charged_trna_counts)
	uncharged_trna.countsIs(uncharged_trna_counts)

def initialize_kinetic_trna_charging(condition, bulk_molecules, sim_data, states):
	'''
	Initializes tRNAs involved in the KineticTrnaChargingModel.
	'''
	def get_f_free(trna, condition):
		if trna == 'selC-tRNA[c]':
			return 0.2 # Initialized to 80% charged

		key = f'{trna}__{condition}'
		if key not in sim_data.relation.trna_condition_to_free_fraction:
			key = f'{trna}__basal'
		return sim_data.relation.trna_condition_to_free_fraction[key]

	# Molecules
	rna_data = sim_data.process.transcription.rna_data
	free_trnas = rna_data['id'][rna_data['is_tRNA']]
	charged_trnas = sim_data.process.transcription.charged_trna_names

	# Fraction free
	f_free = [get_f_free(trna, condition) for trna in free_trnas]

	# Views
	free_view = bulk_molecules.countsView(free_trnas)
	charged_view = bulk_molecules.countsView(charged_trnas)

	# Distribute tRNA counts
	n_total = free_view.counts() + charged_view.counts()
	n_free = np.round(f_free * n_total)
	n_charged = n_total - n_free

	# Update bulk molecules
	free_view.countsIs(n_free)
	charged_view.countsIs(n_charged)

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
