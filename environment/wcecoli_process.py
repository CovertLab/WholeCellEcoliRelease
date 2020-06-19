from __future__ import absolute_import, division, print_function

import os
import copy
import shutil

from vivarium.core.process import Process
from vivarium.core.composition import process_in_experiment
from vivarium.library.dict_utils import deep_merge


from models.ecoli.sim.simulation import ecoli_simulation
from wholecell.utils import constants
import wholecell.utils.filepath as fp
from models.ecoli.sim.variants import apply_variant


DEFAULT_EMITTER = {
	'type': 'database',
	'host': 'localhost:27017',
	'database': 'simulations'}


def initialize_ecoli(config):
	'''
	Args:
		config (dict): options for initializing a simulation
	Returns:
		simulation (CellSimulation): The actual simulation which will perform the calculations.
	'''

	config = copy.deepcopy(config)

	outer_id = config['outer_id']
	sim_out_path = config['sim_out_path']
	variant_type = config['variant_type']
	variant_index = config['variant_index']
	cell_id = config['cell_id']

	sim_data_fit = os.path.join(
		sim_out_path, 'manual', 'kb', constants.SERIALIZED_SIM_DATA_FILENAME)
	variant_sim_data_directory = fp.makedirs(
		sim_out_path, 'agent', outer_id, 'kb', cell_id)
	variant_sim_data_modified_file = os.path.join(
		variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

	# copy sim_data into the experiment directory to support analysis
	shutil.copy(sim_data_fit, variant_sim_data_modified_file)
	fp.verify_file_exists(sim_data_fit, 'Run runParca?')

	# Apply the variant to transform simData.cPickle
	info, sim_data_modified = apply_variant.apply_variant(
		sim_data_file=sim_data_fit,
		variant_type=variant_type,
		variant_index=variant_index)

	config['initialTime'] = config.pop('time') if config.get('time') else 0
	config['simData'] = sim_data_modified
	return ecoli_simulation(**config)

def ecoli_boot_config(agent_config):
	'''
	Instantiates an initial or daughter EcoliSimulation, passes it to a new	`Inner` agent.
	Makes a simOut directory for the simulation in an embedded format:

		out/manual/experiment_id/cohort_id/generation_id/cell_id/simOut

	`agent_config` fields:
		* generation (optional, the cell generation number)
		* outer_id (id of outer environmental agent -- the experiment)
		* working_dir (optional, wcEcoli path containing the sim path out/manual/)
		* files (optional) list of data files:
			files[0] -- inherited_state_path to make a daughter cell
		* start_time (optional)
		* variant_type (optional)
		* variant_index (optional)
		* seed (optional)

	Returns:
		options (dict): simulation arguments for ecoli
	'''
	generation = agent_config.get('generation', 0)
	working_dir = agent_config.get('working_dir', os.getcwd())
	outer_id = agent_config.get('outer_id', 'lattice_000000')
	start_time = agent_config.get('start_time', 0)
	files = agent_config.get('files', [])
	inherited_state_path = files[0] if files else None
	variant_type = agent_config.get('variant_type', 'wildtype')
	variant_index = agent_config.get('variant_index', 0)
	seed = agent_config.get('seed', 0)
	volume = agent_config.get('volume', 1.0)
	cell_id = agent_config.get('cell_id')
	tagged_molecules = agent_config.get('tagged_molecules', ['CDPDIGLYSYN-MONOMER[i]']) # default tag cdsA protein
	emitter_config = agent_config.get('emitter_config', DEFAULT_EMITTER)
	to_report = agent_config.get(
		'to_report',
		{
			'bulk_molecules': [],
			'unique_molecules': [],
			'listeners': [],
		}
	)

	# initialize state
	state = {
		'volume': volume,
		'environment_change': {}}
	agent_config['declare'] = state

	# TODO -- get cohort id (initial cell_id) from lineage trace
	# TODO -- change analysis scripts to allow the cell_id to be used. analysis scripts require starting with 0
	cohort_id = '%06d' % 0
	generation_id = 'generation_%06d' % generation

	# add experiment_id and simulation_id to emitter config (as in lattice_compartment)
	emitter_config['experiment_id'] = outer_id
	emitter_config['simulation_id'] = cell_id

	# make options for boot config
	sim_out_path = fp.makedirs(working_dir, 'out')

	output_dir = os.path.join(sim_out_path, 'agent', outer_id, cohort_id, generation_id, cell_id, 'simOut')
	metadata_dir = fp.makedirs(sim_out_path, 'agent', 'metadata')
	metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)

	options = {
		"sim_out_path":           sim_out_path,
		"variant_type":           variant_type,
		"variant_index":          variant_index,
		"outer_id":               outer_id,
		"lengthSec":              3*60*60,
		"outputDir":              output_dir,
		"initialTime":            start_time,
		"inheritedStatePath":     inherited_state_path,
		"logToDisk":              True,
		"overwriteExistingFiles": True,
		"seed":                   seed,
		"timeStepSafetyFraction": 1.3,
		"maxTimeStep":            0.9,
		"updateTimeStepFreq":     5,
		"logToShell":             True,
		"logToDiskEvery":         1,
		"massDistribution":       True,
		"growthRateNoise":        False,
		"dPeriodDivision":        False,
		"translationSupply":      True,
		"tagged_molecules":       tagged_molecules,
		"to_report":              to_report,
		"emitter_config":         emitter_config,
		"cell_id":                cell_id,
	}

	# Write a metadata file to aid analysis plots.
	metadata = {
		"git_hash":           fp.run_cmdline("git rev-parse HEAD"),
		"git_branch":         fp.run_cmdline("git symbolic-ref --short HEAD"),
		"description":        "an Ecoli Cell Agent",
		"time":               fp.timestamp(),
		"total_gens":         0,  # not known in advance for multi-scale sims
		"analysis_type":      None,
		"variant":            variant_type,
		"mass_distribution":  options['massDistribution'],
		"growth_rate_noise":  options['growthRateNoise'],
		"d_period_division":  options['dPeriodDivision'],
		"translation_supply": options['translationSupply']}
	fp.write_json_file(metadata_path, metadata)

	return options


class wcEcoliAgent(Process):
	defaults = {
		'agent_id': 'X',
		'agent_config': {},
		'time_step': 5.0,  # 60 for big experiments
		'unique_molecules_to_report': [],
		'bulk_molecules_to_report': [],
		'listeners_to_report': [],
	}

	def __init__(self, initial_parameters=None):
		if initial_parameters is None:
			initial_parameters = {}

		parameters = copy.deepcopy(self.defaults)
		deep_merge(parameters, initial_parameters)

		self.agent_id = parameters['agent_id']
		self.agent_config = parameters['agent_config']
		self.agent_config['cell_id'] = self.agent_id
		self.agent_config['working_dir'] = '../wcEcoli'

		self.ecoli_config = ecoli_boot_config(self.agent_config)
		self.ecoli_simulation = initialize_ecoli(self.ecoli_config)
		self.sim_data = self.ecoli_simulation._simData

		environment = self.ecoli_simulation.external_states['Environment']
		media_molecules = set()
		for molecules in environment.saved_media.values():
			media_molecules |= set(molecules)
		self.all_exchange_molecules = list(media_molecules)

		ports_schema = self._ports_schema_from_params(
			parameters, self.all_exchange_molecules
		)
		ports = {
			port: variables_dict.keys()
			for port, variables_dict in ports_schema.items()
		}

		super(wcEcoliAgent, self).__init__(ports, parameters)

	def ports_schema(self):
		return self._ports_schema_from_params(
			self.parameters, self.all_exchange_molecules)

	def _ports_schema_from_params(
		self, parameters, all_exchange_molecules
	):
		ports = [
			'bulk_molecules_report',
			'unique_molecules_report',
			'listeners_report',
			'global',
			'external',
			'exchange',
		]

		schema = {
			port: {} for port in ports
		}

		# local_environment
		schema['external'] = {
			molecule: {
				'_default': 0.0,
				'_emit': True,
				'_updater': 'set',
			}
			for molecule in all_exchange_molecules
		}

		# exchange
		schema['exchange'] = {
			molecule: {
				'_default': 0.0,
			}
			for molecule in all_exchange_molecules
		}

		# bulk_molecules_report
		schema['bulk_molecules_report'] = {
			mol: {
				'_default': None,
				'_emit': True,
				'_updater': 'set',
			}
			for mol in parameters['bulk_molecules_to_report']
		}

		# unique_molecules_report
		schema['unique_molecules_report'] = {
			mol: {
				'_default': None,
				'_emit': True,
				'_updater': 'set',
			}
			for mol in parameters['unique_molecules_to_report']
		}

		# listeners_report
		schema['listeners_report'] = {
			(listener, attr): {
				'_default': None,
				'_emit': True,
				'_updater': 'set',
			}
			for listener, attr in parameters['listeners_to_report']
		}

		# global
		schema['global'] = {
			'mass': {
				# TODO: units
				'_default': 0.0,
				'_emit': True,
				'_updater': 'set',
			},
			'volume': {
				'_default': 0.0,
				'_emit': True,
				'_updater': 'set',
			},
			'media_id': {
				'_default': 'minimal',
				'_emit': True,
				'_updater': 'set',
			},
			'division': {
				'_default': [],
				'_updater': 'set',
			}
		}
		return schema

	def next_update(self, timestep, states):
		media_id = states['global']['media_id']
		local_environment = {
			'concentrations': states['external'],
			'media_id': media_id,
		}
		self.ecoli_simulation.external_states['Environment'].set_local_environment(
			local_environment)
		self.ecoli_simulation.run_for(timestep)
		update = self.ecoli_simulation.generate_inner_update()
		return {
			'global': {
				'volume': update['volume'],
				'division': update['division'],
			},
			'exchange': update['exchange'],
			'unique_molecules_report': update['unique_molecules_report'],
			'bulk_molecules_report': update['bulk_molecules_report'],
			'listeners_report': update['listeners_report'],
		}
