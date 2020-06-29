from __future__ import absolute_import, division, print_function

import os
import copy
import shutil
import sys

from vivarium.core.process import Process
from vivarium.library.dict_utils import deep_merge
from vivarium.library.units import units

from models.ecoli.sim.simulation import ecoli_simulation
from models.ecoli.sim.variants import apply_variant
from wholecell.utils import constants
import wholecell.utils.filepath as fp


def initialize_ecoli(config):
	'''Create a Simulation object

	Instantiates an initial or daughter EcoliSimulation, and returns it.
	Copies the simData file to a new directory for the instantiated
	agent at: ``out/agent/outer_id/kb/agent_id``.

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
	_, sim_data_modified = apply_variant.apply_variant(
		sim_data_file=sim_data_fit,
		variant_type=variant_type,
		variant_index=variant_index)

	config['initialTime'] = config.pop('time') if config.get('time') else 0
	config['simData'] = sim_data_modified
	return ecoli_simulation(**config)


def ecoli_boot_config(agent_config):
	'''Generate a configuration for :py:func:`intialize_ecoli`

	Makes a simOut directory for the simulation at:
	``out/manual/experiment_id/cohort_id/generation_id/cell_id/simOut``

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
	* to_report (optional): A dictionary specifying lists of bulk
	  molecules, unique molecules, and listener attributes to include in
	  updates. The ``bulk_molecules`` key should be mapped to a list of
	  bulk molecule names. The ``unique_molecules`` key should similarly
	  be mapped to a list of unique molecule names. The ``listeners``
	  key should be mapped to a list of tuples, each of which contains
	  as the first element the name of the listener and as the second
	  key the attribute of that listener with the value to report.

	Returns:
		options (dict): simulation arguments for ecoli
	'''
	generation = agent_config.get('generation', 0)
	working_dir = agent_config.get('working_dir', os.getcwd())
	outer_id = agent_config.get('outer_id', 'lattice_000000')
	length_sec = agent_config.get('lengthSec', 3*60*60)
	start_time = agent_config.get('start_time', 0)
	files = agent_config.get('files', [])
	inherited_state_path = files[0] if files else None
	variant_type = agent_config.get('variant_type', 'wildtype')
	variant_index = agent_config.get('variant_index', 0)
	seed = agent_config.get('seed', 0)
	volume = agent_config.get('volume', 1.0)
	cell_id = agent_config.get('cell_id')
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
		'environment_change': {},
	}
	agent_config['declare'] = state

	# TODO -- get cohort id (initial cell_id) from lineage trace
	cohort_id = '%06d' % 0
	generation_id = 'generation_%06d' % generation

	# make options for boot config
	sim_out_path = fp.makedirs(working_dir, 'out')

	output_dir = os.path.join(
		sim_out_path, 'agent', outer_id, cohort_id, generation_id, cell_id,
		'simOut'
	)
	metadata_dir = fp.makedirs(sim_out_path, 'agent', 'metadata')
	metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)

	options = {
		"sim_out_path":           sim_out_path,
		"variant_type":           variant_type,
		"variant_index":          variant_index,
		"outer_id":               outer_id,
		"lengthSec":              length_sec,
		"outputDir":              output_dir,
		"initialTime":            start_time,
		"inheritedStatePath":     inherited_state_path,
		"logToDisk":              False,  # We're using Vivarium instead
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
		"to_report":              to_report,
		"cell_id":                cell_id,
	}

	# Write a metadata file to aid analysis plots.
	metadata = {
		"git_hash":           fp.run_cmdline("git rev-parse HEAD"),
		"git_branch":         fp.run_cmdline("git symbolic-ref --short HEAD"),
		"description":        "an Ecoli Cell Agent",
		"time":               fp.timestamp(),
		"python":             sys.version.splitlines()[0],
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
		'agent_config': {
			'to_report': {
				'listeners': [
					('Mass', 'cellMass'),
					('Mass', 'cellDensity'),
					('Mass', 'volume'),
				],
				'unique_molecules': [],
				'bulk_molecules': [],
			},
		},
		'time_step': 5.0,  # 60 for big experiments
	}
	# Must match units used by wcEcoli
	mass_units = units.fg

	def __init__(self, initial_parameters=None):
		'''Process that internally runs a wcEcoli simulation

		This process wraps a wcEcoli simulation, passing it information
		about the environment and returning its updates. It also reports
		parts of the internal state of the wcEcoli simulation.

		.. WARNING:: If you specify listeners to report, you MUST
			include the following tuples:

			* ``('Mass', 'cellMass')``
			* ``('Mass', 'cellDensity')``
			* ``('Mass', 'volume')``

			These listeners are used by this process to update the
			``globals`` port.

		Ports:

		* **bulk_molecules_report**: Holds a variable for each bulk
		  molecule to be reported from wcEcoli. This variable will be
		  updated to hold the current count of the molecule from
		  wcEcoli. Should be mapped to the agent's boundary.
		* **unique_molecules_report**: Holds a variable for each unique
		  molecule to be reported from wcEcoli. This variable will be
		  updated to hold the current count of the molecule from
		  wcEcoli. Should be mapped to the agent's boundary.
		* **listeners_report**: Holds a variable for each listener
		  attribute to be reported from wcEcoli. For a listener
		  ``listen`` with attribute ``attr``, the value of
		  ``listen.attr`` will be reported as the value of the variable
		  ``listen-attr``. Should be mapped to the agent's boundary.
		* **global**: Should be mapped to the agent's boundary.
		* **external**: Should be updated by other processes to contain
		  the local molecule concentrations around this agent.
		* **exchange**: Will reflect the flux of molecules into or out
		  of the cell.

		Arguments:
			initial_parameters (dict): Accepts the following
				configuration keys:

				* **agent_id** (:py:class:`str`): Uniquely identifies
				  this agent within the Vivarium model.
				* **agent_config** (:py:class:`dict`): Configuration
				  options that will be passed to
				  :py:func:`ecoli_boot_config` when generating the
				  wcEcoli object.
		'''
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
		self.sim_data = self.ecoli_simulation.get_sim_data()

		environment = self.ecoli_simulation.external_states['Environment']
		media_molecules = set()
		for molecules in environment.saved_media.values():
			media_molecules |= set(molecules)
		self.all_exchange_molecules = list(media_molecules)

		super(wcEcoliAgent, self).__init__(parameters)

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
		if (
			parameters['agent_config']['to_report']['bulk_molecules']
		) == []:
			schema['bulk_molecules_report'] = {
				'_default': {},
				'_emit': True,
				'_updater': 'set',
			}
		else:
			schema['bulk_molecules_report'] = {
				mol: {
					'_default': 0.0,
					'_emit': True,
					'_updater': 'set',
				}
				for mol
				in parameters['agent_config']['to_report']['bulk_molecules']
			}

		# unique_molecules_report
		if (
			parameters['agent_config']['to_report']['unique_molecules']
		) == []:
			schema['unique_molecules_report'] = {
				'_default': {},
				'_emit': True,
				'_updater': 'set',
			}
		else:
			schema['unique_molecules_report'] = {
				mol: {
					'_default': 0.0,
					'_emit': True,
					'_updater': 'set',
				}
				for mol
				in parameters['agent_config']['to_report']['unique_molecules']
			}

		# listeners_report
		if (
			parameters['agent_config']['to_report']['listeners']
		) == []:
			schema['listeners_report'] = {
				'_default': {},
				'_emit': True,
				'_updater': 'set',
			}
		else:
			schema['listeners_report'] = {
				'-'.join((listener, attr)): {
					'_default': 0.0,
					'_emit': True,
					'_updater': 'set',
				}
				for listener, attr
				in parameters['agent_config']['to_report']['listeners']
			}

		# global
		schema['global'] = {
			'mass': {
				'_default': 0.0 * self.mass_units,
				'_emit': True,
				'_updater': 'set',
			},
			'volume': {
				'_default': 0.0,
				'_emit': True,
				'_updater': 'set',
			},
			'density': {
				'_default': 1.0,
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
		listeners_report = update['listeners_report']
		return {
			'global': {
				'volume': (
					listeners_report[('Mass', 'volume')]
				),
				'mass': (
					listeners_report[('Mass', 'cellMass')]
					* self.mass_units
				),
				'density': (
					listeners_report[('Mass', 'cellDensity')]
				),
				'division': update['division'],
			},
			'exchange': update['exchange'],
			'unique_molecules_report': update['unique_molecules_report'],
			'bulk_molecules_report': update['bulk_molecules_report'],
			'listeners_report': {
				'-'.join(key): value
				for key, value in listeners_report.items()
			},
		}
