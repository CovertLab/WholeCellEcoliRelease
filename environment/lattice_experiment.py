from __future__ import absolute_import, division, print_function

import json
from urllib import quote_plus

from vivarium.core.composition import (
	make_agents,
	make_experiment_from_compartment_dicts,
	simulate_experiment,
	plot_agents_multigen,
)
from vivarium.compartments.lattice import Lattice
from vivarium.plots.multibody_physics import plot_snapshots

from wcecoli_process import wcEcoliAgent
from wcecoli_compartment import WcEcoliCell


MEDIA_IDS = ["minimal", "minimal_minus_oxygen",
		"minimal_plus_amino_acids"]
MEDIA_ID = "minimal"
N_WCECOLI_AGENTS = 1  # Works at 50
BOUNDS = (50, 50)
N_BINS = (20, 20)
OUT_DIR = 'out/agent'
SAVE_OUTPUT = False
SIMULATION_TIME = 20
DEFAULT_EMITTER_CONFIG = {
	'type': 'database',
	'host': 'localhost:27017',
	'database': 'simulations',
}
SECRETS_PATH = 'environment/secrets.json'
ENVIRONMENT_CONFIG_KEYS = ['bounds', 'size', 'agents']


def get_atlas_database_emitter_config(
	username, password, cluster_subdomain, database
):
	username = quote_plus(username)
	password = quote_plus(password)
	database = quote_plus(database)

	uri = (
		"mongodb+srv://{}:{}@{}.mongodb.net/"
		+ "?retryWrites=true&w=majority"
	).format(username, password, cluster_subdomain)
	return {
		'type': 'database',
		'host': uri,
		'database': database,
	}


def emit_environment_config(environment_config, emitter):
	config = {
        'bounds': environment_config['multibody']['bounds'],
		'type': 'environment_config',
	}
	emitter.emit({
		'data': config,
		'table': 'configuration',
	})


def main():
	with open(SECRETS_PATH, 'r') as f:
		secrets = json.load(f)
	emitter_config = get_atlas_database_emitter_config(
		**secrets['database'])
	# To use a timeseries emitter, uncomment the line below
	# emitter_config = {'type': 'timeseries'}
	agent = wcEcoliAgent({})
	external_states = agent.ecoli_simulation.external_states
	# Assert agent has media_id MEDIA_ID
	recipe = external_states['Environment'].saved_media[MEDIA_ID]

	compartment = WcEcoliCell()
	agent_ids = ['wcecoli_{}'.format(i) for i in range(N_WCECOLI_AGENTS)]
	agents_dict = make_agents(agent_ids, compartment)

	environment_config = {
		'multibody': {
			'bounds': BOUNDS,
			'size': BOUNDS,
			'agents': agents_dict,
		},
		'diffusion': {
			'bounds': BOUNDS,
			'n_bins': N_BINS,
			'molecules': recipe.keys(),
            'depth': 10000.0,  # Deep to avoid depleting local molecules
			'diffusion': 5,  # 10x faster than the default 5e-1
			'gradient': {
				'type': 'linear',
				'molecules': {
					molecule: {
						'center': (0, 0),
						'base': concentration,
						'slope': 0,
					}
					for molecule, concentration in recipe.items()
				},
			},
		},
	}
	environment = Lattice(environment_config)
	environment_dict = environment.generate()

	initial_state = {}
	experiment = make_experiment_from_compartment_dicts(
		environment_dict, agents_dict, emitter_config, initial_state)
	print('Experiment ID:', experiment.experiment_id)
	emit_environment_config(environment_config, experiment.emitter)
	settings = {
		'timestep': 1.0,
		'total_time': SIMULATION_TIME,
		'return_raw_data': True,
	}
	data = simulate_experiment(experiment, settings)

	if SAVE_OUTPUT:
		# Plot snapshots
		agents = {
			time: time_data['agents'] for time, time_data in data.items()}
		fields = {
			time: time_data['fields'] for time, time_data in data.items()}
		snapshots_data = {
			'agents': agents,
			'fields': fields,
			'config': environment_config['multibody'],
		}
		plot_config = {
			'out_dir': OUT_DIR,
			'filename': 'snapshot',
		}
		plot_snapshots(snapshots_data, plot_config)

		# Plot Emitted Values
		plot_settings = {
			'agents_key': 'agents',
		}
		plot_agents_multigen(data, plot_settings, OUT_DIR, 'simulation')


if __name__ == '__main__':
	main()
