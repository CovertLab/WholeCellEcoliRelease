import numpy as np

from models.ecoli.processes.polypeptide_elongation import SteadyStateElongationModel
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
	init.initializeUniqueMoleculesFromBulk(bulkMolCntr, uniqueMolCntr, sim_data,
		randomState, sim._superhelical_density, sim._trna_attenuation)

	# Must be called after unique and bulk molecules are initialized to get
	# concentrations for ribosomes, tRNA, synthetases etc from cell volume
	if sim._trna_charging:
		elongation_model = SteadyStateElongationModel(sim_data, sim.processes['PolypeptideElongation'])
		initialize_trna_charging(sim_data, sim.internal_states, elongation_model.calculate_trna_charging)

	# Adjust small molecule concentrations again after other mass adjustments
	# for more stable metabolism solution at beginning of sims
	init.set_small_molecule_counts(bulkMolCntr, sim_data, media_id, import_molecules,
		massCoeff, cell_mass=init.calculate_cell_mass(sim.internal_states))

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
	fraction_charged, _ = calc_charging(synthetase_conc, uncharged_trna_conc, charged_trna_conc, aa_conc, ribosome_conc, f)

	# Update counts of tRNA to match charging
	total_trna_counts = uncharged_trna.counts() + charged_trna.counts()
	charged_trna_counts = np.round(total_trna_counts * np.dot(fraction_charged, aa_from_trna))
	uncharged_trna_counts = total_trna_counts - charged_trna_counts
	charged_trna.countsIs(charged_trna_counts)
	uncharged_trna.countsIs(uncharged_trna_counts)

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
