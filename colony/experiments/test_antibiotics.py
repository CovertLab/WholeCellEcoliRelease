'''Test script for the Antibiotics Experiment of an E. coli colony
'''

from __future__ import absolute_import, division, print_function

import argparse
import os

from vivarium.core.composition import (
	load_timeseries,
	assert_timeseries_close
)
from vivarium.library.timeseries import (
	process_path_timeseries_for_csv,
	save_flat_timeseries,
)

from colony.experiments.antibiotics import simulate
from colony.constants import REFERENCE_DATA_PATH
from colony.constants import OUT_DIR


#: Location of the reference data
REFERENCE_DATA_FILE = os.path.join(
	REFERENCE_DATA_PATH, 'antibiotics.csv')
#: File to write simulation data to
OUTPUT_DATA_PATH = os.path.join(
	OUT_DIR, 'antibiotics_test_output.csv')
#: Length of simulation to run for both generation and checking
DURATION = 1300
#: Number of cells with which to initialize simulation
NUM_CELLS = 1
PULSE_CONCENTRATION = 0
ADD_AA = True
ANAEROBIC = False
ANTIBIOTIC_THRESHOLD = 0.02
UPDATE_FIELDS = False


def run_simulation():
	# type: () -> dict
	'''Run a simulation

	Returns:
		dict: A processed path timeseries. Processing is done with
		:py:func:`process_path_timeseries_for_csv`
	'''
	emitter_config = {
		'type': 'timeseries'
	}
	emitter = simulate(
		emitter_config, DURATION, NUM_CELLS, PULSE_CONCENTRATION,
		ADD_AA, ANAEROBIC, ANTIBIOTIC_THRESHOLD, UPDATE_FIELDS,
	)
	path_ts = emitter.get_path_timeseries()
	processed_ts = process_path_timeseries_for_csv(path_ts)
	return processed_ts


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'--generate', '-g',
		action='store_true',
		default=False,
		help='Save reference data to {}.'.format(REFERENCE_DATA_FILE),
	)
	parser.add_argument(
		'--check', '-c',
		action='store_true',
		default=False,
		help='Check that model behavior matches {}.'.format(
            REFERENCE_DATA_FILE),
	)
	args = parser.parse_args()
	assert args.check != args.generate
	processed_ts = run_simulation()
	save_flat_timeseries(
		processed_ts,
		out_dir=os.path.dirname(OUTPUT_DATA_PATH),
		filename=os.path.basename(OUTPUT_DATA_PATH),
	)
	if args.check:
		reference_ts = load_timeseries(REFERENCE_DATA_FILE)
		# We load simulation data from the file so that all
		# transformations by CSV reader/writer are applied
		simulation_ts = load_timeseries(OUTPUT_DATA_PATH)
		assert_timeseries_close(
			simulation_ts, reference_ts,
			tolerances={
				'agents,wcecoli_0,boundary,bulk_molecules_report,BCCP-MONOMER[c]':
					900,
				'agents,wcecoli_0,boundary,bulk_molecules_report,G7602-MONOMER[c]':
					99,
				'agents,wcecoli_0,boundary,bulk_molecules_report,G7763-MONOMER[c]':
					9,
				'agents,wcecoli_0,boundary,bulk_molecules_report,EG11162-MONOMER[c]':
					50,
				'agents,wcecoli_0,boundary,bulk_molecules_report,EG11256-MONOMER[c]':
					500,
				'agents,wcecoli_0,boundary,dry_mass': 9,
				'agents,wcecoli_0,boundary,mass': 9,
				'agents,wcecoli_0,boundary,bulk_molecules_report,PD03585[c]':
					9,
			},
		)
	if args.generate:
		save_flat_timeseries(
			processed_ts,
			os.path.dirname(OUTPUT_DATA_PATH),
			os.path.basename(OUTPUT_DATA_PATH),
		)


if __name__ == '__main__':
	main()
