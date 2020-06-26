'''Simulate a colony of E. coli

Each E. coli cell is modeled with wcEcoli. They are placed into a shared
environment using Vivarium.
'''

from __future__ import absolute_import, division, print_function

import argparse
import csv
import json

from vivarium.core.composition import (
	make_agents,
	make_experiment_from_compartment_dicts,
	simulate_experiment,
)
from vivarium.compartments.lattice import Lattice
from vivarium.core.emitter import (
	get_atlas_database_emitter_config,
	emit_environment_config,
	SECRETS_PATH,
)

from environment.wcecoli_process import wcEcoliAgent
from environment.wcecoli_compartment import WcEcoliCell


MEDIA_IDS = ["minimal", "minimal_minus_oxygen",
		"minimal_plus_amino_acids"]
MEDIA_ID = "minimal"
N_WCECOLI_AGENTS = 1  # Works at 50
BOUNDS = (50, 50)
N_BINS = (20, 20)
DEFAULT_SIMULATION_TIME = 60 * 60 * 1.5  # 1.5 hr
TAGGED_MOLECULES_PATH = 'environment/tagged_molecules.csv'


def simulate(emitter_config, simulation_time):
	'''Run the simulation

	Arguments:
		emitter_config: See the Vivarium docs for details. This
			specifies how the simulation data is reported, e.g. to a
			database.
		simulation_time: Seconds of time to simulate.

	Returns:
		vivarium.core.emitter.Emitter: An emitter from which the
		simulation data can be extracted. If you choose to emit to a
		database, you can ignore the returned emitter as the data is
		sent automatically.
	'''
	agent = wcEcoliAgent({})
	external_states = agent.ecoli_simulation.external_states
	# TODO: Assert agent has media_id MEDIA_ID
	recipe = external_states['Environment'].saved_media[MEDIA_ID]

	with open(TAGGED_MOLECULES_PATH, 'r') as f:
		reader = csv.reader(f)
		tagged_molecules = [
			molecule for _, molecule in reader
		]
	compartment = WcEcoliCell()
	agent_ids = ['wcecoli_{}'.format(i) for i in range(N_WCECOLI_AGENTS)]
	process_config = {
		'agent_config': {
			'to_report': {
				'bulk_molecules': tagged_molecules,
			},
		},
	}
	agents_dict = make_agents(agent_ids, compartment, process_config)

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
		'total_time': simulation_time,
		'return_raw_data': True,
	}
	simulate_experiment(experiment, settings)
	return experiment.emitter


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'--atlas', '-a',
		action='store_true',
		default=False,
		help=(
			'Read data from an mongoDB Atlas instead of a local mongoDB. '
			'Credentials, cluster subdomain, and database name should be '
			'specified in {}.'.format(SECRETS_PATH)
		)
	)
	parser.add_argument(
		'--port', '-p',
		default=27017,
		type=int,
		help='Port at which to access local mongoDB instance.',
	)
	parser.add_argument(
		'--host', '-o',
		default='localhost',
		type=str,
		help='Host at which to access local mongoDB instance.',
	)
	parser.add_argument(
		'--database_name', '-d',
		default='simulations',
		type=str,
		help=(
			'Name of database on local mongoDB instance to read from.'
		)
	)
	parser.add_argument(
		'--simulation_time', '-s',
		default=DEFAULT_SIMULATION_TIME,
		type=int,
		help='Number of seconds to simulate.',
	)
	args = parser.parse_args()
	if args.atlas:
		with open(SECRETS_PATH, 'r') as f:
			secrets = json.load(f)
		emitter_config = get_atlas_database_emitter_config(
			**secrets['database'])
	else:
		emitter_config = {
			'type': 'database',
			'host': '{}:{}'.format(args.host, args.port),
			'database': args.database_name,
		}
	_ = simulate(emitter_config, args.simulation_time)


if __name__ == '__main__':
	main()
