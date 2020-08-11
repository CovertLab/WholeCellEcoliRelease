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
from vivarium.core.composition import (
	agent_environment_experiment,
	simulate_experiment,
)
from vivarium.compartments.lattice import Lattice
from vivarium.core.emitter import (
	get_atlas_database_emitter_config,
	emit_environment_config,
	SECRETS_PATH,
)
from vivarium.processes.diffusion_field import make_gradient

from wholecell.io import tsv
from environment.wcecoli_process import wcEcoliAgent
from environment.wcecoli_compartment import (
	ANTIBIOTIC_KEY,
	INITIAL_EXTERNAL_ANTIBIOTIC,
	WcEcoliCell,
)


MINIMAL_MEDIA_ID = 'minimal'
AA_MEDIA_ID = 'minimal_plus_amino_acids'
#: From simData.ordered_conditions
CONDITION_VARIANTS = [
	'basal', 'no_oxygen', 'with_aa', 'acetate', 'succinate']
BOUNDS = (50, 50)
N_BINS = (20, 20)
TAGGED_MOLECULES_PATH = os.path.join(
	os.path.dirname(__file__), 'tagged_molecules.csv')
NUM_EMISSIONS = 100


def get_timeline(n_bins, size, pulses, end_time):
	'''Get a timeline for antibiotic pulses.

	Arguments:
		n_bins (list): Number of bins in x and y directions.
		size (list): Size of environment in x and y directions.
		pulses (list): List of tuples, each of which describes a pulse.
			Each tuple has the form (start_time, duration,
			concentration).
		end_time (int): The length of the experiment
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


class TestGetTimeline:

	def assert_timelines_equal(self, timeline, expected_timeline):
		assert len(timeline) == len(expected_timeline)
		for (actual_time, actual), (expected_time, expected) in zip(
			timeline, expected_timeline
		):
			assert actual_time == expected_time
			if expected == {}:
				assert actual == expected
			else:
				assert (
					actual[('fields', ANTIBIOTIC_KEY)].tolist()
					== expected[('fields', ANTIBIOTIC_KEY)].tolist()
				)
			assert len(actual) == len(expected)

	def test_one_pulse_one_bin(self):
		timeline = get_timeline(
			(1, 1), (1, 1), [(1, 1, 1)], 5)
		expected_timeline = [
			(1, {('fields', ANTIBIOTIC_KEY): np.ones((1, 1))}),
			(2, {('fields', ANTIBIOTIC_KEY): np.zeros((1, 1))}),
			(5, {}),
		]
		assert timeline == expected_timeline

	def test_selects_correct_bin(self):
		timeline = get_timeline(
			(2, 2), (1, 1), [(1, 1, 1)], 5)
		expected_timeline = [
			(1, {('fields', ANTIBIOTIC_KEY): np.ones((2, 2))}),
			(2, {('fields', ANTIBIOTIC_KEY): np.zeros((2, 2))}),
			(5, {}),
		]
		self.assert_timelines_equal(timeline, expected_timeline)

	def test_two_pulses(self):
		timeline = get_timeline(
			(1, 1),
			(1, 1),
			[(1, 1, 1), (4, 2, 2)],
			10,
		)
		expected_timeline = [
			(1, {('fields', ANTIBIOTIC_KEY): np.ones((1, 1))}),
			(2, {('fields', ANTIBIOTIC_KEY): np.zeros((1, 1))}),
			(4, {('fields', ANTIBIOTIC_KEY): np.full((1, 1), 2)}),
			(6, {('fields', ANTIBIOTIC_KEY): np.zeros((1, 1))}),
			(10, {}),
		]
		self.assert_timelines_equal(timeline, expected_timeline)


def simulate(
	emitter_config, simulation_time, num_cells, pulse_concentration,
	add_aa, antibiotic_threshold,
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
		antibiotic_threshold: The maximum internal concentration of
			antibiotic cells can survive.

	Returns:
		vivarium.core.emitter.Emitter: An emitter from which the
		simulation data can be extracted. If you choose to emit to a
		database, you can ignore the returned emitter as the data is
		sent automatically.
	'''
	agent = wcEcoliAgent({})
	external_states = agent.ecoli_simulation.external_states
	media_id = AA_MEDIA_ID if add_aa else MINIMAL_MEDIA_ID
	variant_name = 'with_aa' if add_aa else 'basal'
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
		simulation_time * 0.25,
		pulse_concentration,
	)
	timeline_config = {
		'timeline': get_timeline(
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
		'update_fields': False,
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
		'type': WcEcoliCell,
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
		'--antibiotic_threshold',
		type=float,
		default=0.86,
		help='Internal antibiotic concentration past which cells die.',
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
		args.antibiotic_threshold,
	)


if __name__ == '__main__':
	main()
