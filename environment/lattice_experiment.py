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
from vivarium.library.lattice_utils import get_bin_site

from wholecell.io import tsv
from environment.wcecoli_process import wcEcoliAgent
from environment.wcecoli_compartment import (
	ANTIBIOTIC_KEY,
	INITIAL_EXTERNAL_ANTIBIOTIC,
	WcEcoliCell,
)


MEDIA_ID = "minimal"
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
			Each tuple has the form
			(start_time, duration, center, concentration) where center
			is a 2-tuple that specifies coordinates for the pulse.
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
	for start_time, duration, center, concentration in pulses:
		center_row, center_col = get_bin_site(center, n_bins, size)
		start_field = get_empty_field()
		start_field[center_row, center_col] = concentration
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
			(1, 1), (1, 1), [(1, 1, (0.5, 0.5), 1)], 5)
		expected_timeline = [
			(1, {('fields', ANTIBIOTIC_KEY): np.ones((1, 1))}),
			(2, {('fields', ANTIBIOTIC_KEY): np.zeros((1, 1))}),
			(5, {}),
		]
		assert timeline == expected_timeline

	def test_selects_correct_bin(self):
		timeline = get_timeline(
			(2, 2), (1, 1), [(1, 1, (0.75, 0.25), 1)], 5)
		expected_timeline = [
			(1, {('fields', ANTIBIOTIC_KEY):
				np.array([[0, 0], [1, 0]])}),
			(2, {('fields', ANTIBIOTIC_KEY): np.zeros((2, 2))}),
			(5, {}),
		]
		self.assert_timelines_equal(timeline, expected_timeline)

	def test_two_pulses(self):
		timeline = get_timeline(
			(1, 1),
			(1, 1),
			[(1, 1, (0.5, 0.5), 1), (4, 2, (0.75, 0.25), 2)],
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


def simulate(emitter_config, simulation_time, num_cells, length_sec=None):
	'''Run the simulation

	Arguments:
		emitter_config: See the Vivarium docs for details. This
			specifies how the simulation data is reported, e.g. to a
			database.
		simulation_time: Seconds of time to simulate.
		num_cells: Number of cells to initialize simulation with.
		length_sec: Seconds before we force division to occur. This
			should be left at its default for simulations. It is useful
			to set it to some small value (e.g. 20 seconds) for testing
			or debugging.

	Returns:
		vivarium.core.emitter.Emitter: An emitter from which the
		simulation data can be extracted. If you choose to emit to a
		database, you can ignore the returned emitter as the data is
		sent automatically.
	'''
	agent = wcEcoliAgent({})
	external_states = agent.ecoli_simulation.external_states
	recipe = external_states['Environment'].saved_media[MEDIA_ID]
	recipe[ANTIBIOTIC_KEY] = INITIAL_EXTERNAL_ANTIBIOTIC

	with io.open(TAGGED_MOLECULES_PATH, 'rb') as f:
		reader = tsv.reader(f, delimiter=',')
		tagged_molecules = [
			molecule for _, _, molecule in reader
		]
	antibiotic_pulse = (
		simulation_time * 0.5,
		simulation_time * 0.25,
		[0.5 * bound for bound in BOUNDS],
		1000,
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
		},
		'_parallel': True,
		'update_fields': False,
		'timeline': timeline_config,
	}
	if length_sec is not None:
		process_config['agent_config']['lengthSec'] = length_sec
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
		'--length_sec', '-l',
		default=-1,
		type=int,
		help='Seconds until cell division is forced.',
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
	length_sec = args.length_sec if args.length_sec >= 0 else None
	_ = simulate(
		emitter_config,
		args.simulation_time,
		args.num_cells,
		length_sec
	)


if __name__ == '__main__':
	main()
