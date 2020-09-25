'''Simulate a colony of E. coli

Each E. coli cell is modeled with wcEcoli. They are placed into a shared
environment using Vivarium.
'''

from __future__ import absolute_import, division, print_function

import argparse
import io
import json
import os

import numpy as np
from typing import List, Tuple
from vivarium.core.composition import (
	agent_environment_experiment,
	simulate_experiment,
)
from vivarium.core.emitter import (
	get_atlas_database_emitter_config,
	emit_environment_config,
	SECRETS_PATH,
)
from vivarium_cell.composites.lattice import Lattice
from vivarium_cell.processes.diffusion_field import make_gradient

from wholecell.io import tsv
from colony.processes.wcecoli import WcEcoli
from colony.compartments.antibiotics import (
	ANTIBIOTIC_KEY,
	INITIAL_EXTERNAL_ANTIBIOTIC,
	AntibioticsCell,
)


MINIMAL_MEDIA_ID = 'minimal'
AA_MEDIA_ID = 'minimal_plus_amino_acids'
ANAEROBIC_MEDIA_ID = 'minimal_minus_oxygen'
#: From simData.ordered_conditions
CONDITION_VARIANTS = [
	'basal', 'no_oxygen', 'with_aa', 'acetate', 'succinate']
BOUNDS = (50, 50)
N_BINS = (20, 20)
TAGGED_MOLECULES_PATH = os.path.join(
	os.path.dirname(__file__), 'antibiotics_tagged_molecules.csv')
NUM_EMISSIONS = 100


def get_antibiotics_timeline(
	n_bins,  # type: Tuple[int, int]
	size,  # type: Tuple[int, int]
	pulses,  # type: List[Tuple[int, int, float]]
	end_time,  # type: int
):
	'''Get a timeline for antibiotic pulses.

	Arguments:
		n_bins: Number of bins in x and y directions.
		size: Size of environment in x and y directions.
		pulses: List of tuples, each of which describes a pulse.
			Each tuple has the form (start_time, duration,
			concentration).
		end_time: The length of the experiment
	Returns:
		list: A timeline that implements the described pulses.
	'''
	timeline = []
	get_empty_field = lambda: make_gradient(
		{
			'type': 'uniform',
			'molecules': {
				ANTIBIOTIC_KEY: '0',
			},
		},
		n_bins,
		size,
	)[ANTIBIOTIC_KEY]
	for start_time, duration, concentration in pulses:
		start_field = np.full_like(get_empty_field(), concentration)
		timeline.append((
			start_time,
			{('fields', ANTIBIOTIC_KEY): start_field},
		))
		end_field = get_empty_field()
		timeline.append((
			start_time + duration,
			{('fields', ANTIBIOTIC_KEY): end_field},
		))
	timeline.append((end_time, {}))
	return timeline


def simulate(
	emitter_config, simulation_time, num_cells, pulse_concentration,
	add_aa, anaerobic, antibiotic_threshold, update_fields,
):
	'''Run the simulation

	Arguments:
		emitter_config: See the Vivarium docs for details. This
			specifies how the simulation data is reported, e.g. to a
			database.
		simulation_time: Seconds of time to simulate.
		num_cells: Number of cells to initialize simulation with.
		pulse_concentration: Concentration of antibiotic to provide
			during pulse.
		add_aa: Whether to add amino acids to the media and use the
			with_aa wcEcoli variant.
		anaerobic: Whether to run an aerobic variant in anaerobic media.
		antibiotic_threshold: The maximum internal concentration of
			antibiotic cells can survive.
		update_fields: Whether to let wcEcoli update environmental
			fields. This often needs to be set to False to avoid
			breaking FBA.

	Returns:
		vivarium.core.emitter.Emitter: An emitter from which the
		simulation data can be extracted. If you choose to emit to a
		database, you can ignore the returned emitter as the data is
		sent automatically.
	'''
	agent = WcEcoli({})
	external_states = agent.ecoli_simulation.external_states
	if add_aa and anaerobic:
		raise ValueError(
			'At most one of add_aa and anaerobic may be True.')
	elif add_aa:
		media_id = AA_MEDIA_ID
		variant_name = 'with_aa'
	elif anaerobic:
		media_id = ANAEROBIC_MEDIA_ID
		variant_name = 'no_oxygen'
	else:
		media_id = MINIMAL_MEDIA_ID
		variant_name = 'basal'
	variant_index = CONDITION_VARIANTS.index(variant_name)
	recipe = external_states['Environment'].saved_media[media_id]
	recipe[ANTIBIOTIC_KEY] = INITIAL_EXTERNAL_ANTIBIOTIC

	with io.open(TAGGED_MOLECULES_PATH, 'rb') as f:
		reader = tsv.reader(f, delimiter=',')
		tagged_molecules = [
			molecule for _, _, molecule in reader
		]
	antibiotic_pulse = (
		simulation_time * 0.5,
		simulation_time * 0.5,
		pulse_concentration,
	)
	timeline_config = {
		'timeline': get_antibiotics_timeline(
			N_BINS, BOUNDS, [antibiotic_pulse], simulation_time),
	}
	process_config = {
		'agent_config': {
			'to_report': {
				'bulk_molecules': tagged_molecules,
			},
			'media_id': media_id,
			'variant_type': 'condition',
			'variant_index': variant_index,
		},
		'_parallel': True,
		'update_fields': update_fields,
		'timeline': timeline_config,
		'death': {
			'detectors': {
				'antibiotic': {
					'antibiotic_threshold': antibiotic_threshold,
				},
			},
		},
	}
	agent_ids = ['wcecoli_{}'.format(i) for i in range(num_cells)]

	initial_state = {}
	settings = {
		'emitter': emitter_config,
		'emit_step': max(simulation_time // NUM_EMISSIONS, 1)
	}
	agents_config = {
		'type': AntibioticsCell,
		'ids': agent_ids,
		'config': process_config,
	}
	environment_config = {
		'type': Lattice,
		'config': {
			'multibody': {
				'bounds': BOUNDS,
			},
			'diffusion': {
				'bounds': BOUNDS,
				'n_bins': N_BINS,
				'molecules': list(recipe.keys()),
				'depth': 1000,  # Deep to avoid depleting local molecules
				'diffusion': 5e-1,
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
			'colony_shape_deriver': {
				'bounds': BOUNDS,
				'alpha': 0.2,
			},
		},
	}
	experiment = agent_environment_experiment(
		agents_config,
		environment_config,
		initial_state,
		settings=settings,
	)
	print('Experiment ID:', experiment.experiment_id)
	emit_environment_config(environment_config['config'], experiment.emitter)
	settings = {
		'timestep': 1.0,
		'return_raw_data': True,
		'timeline': timeline_config,
	}
	simulate_experiment(experiment, settings)
	experiment.end()
	print('Experiment ID:', experiment.experiment_id)
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
		default=60 * 60 * 1.5,  # 1.5 hr
		type=int,
		help='Number of seconds to simulate.',
	)
	parser.add_argument(
		'--num_cells', '-n',
		default=1,
		type=int,
		help='Number of cells to create at start of simulation.',
	)
	parser.add_argument(
		'--pulse_concentration',
		default=1,
		type=float,
		help='Antibiotic concentration to provide in pulse.'
	)
	parser.add_argument(
		'--add_aa',
		action='store_true',
		help='Use media with amino acids and the corresponding variant.'
	)
	parser.add_argument(
		'--anaerobic',
		action='store_true',
		help='Use anaerobic media and variant.'
	)
	parser.add_argument(
		'--antibiotic_threshold',
		type=float,
		default=0.86,
		help='Internal antibiotic concentration past which cells die.',
	)
	parser.add_argument(
		'--update_fields',
		action='store_true',
		help='Let wcEcoli update the environment fields.',
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
	_ = simulate(
		emitter_config,
		args.simulation_time,
		args.num_cells,
		args.pulse_concentration,
		args.add_aa,
		args.anaerobic,
		args.antibiotic_threshold,
		args.update_fields,
	)


if __name__ == '__main__':
	main()
