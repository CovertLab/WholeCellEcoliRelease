'''Test script for the Lattice Experiment of an E. coli colony
'''

from __future__ import absolute_import, division, print_function

import argparse
import io
import os

from vivarium.core.composition import (
	load_timeseries,
	assert_timeseries_close
)
from wholecell.io import tsv

from environment.lattice_experiment import simulate


#: Location of the reference data
REFERENCE_DATA_PATH = os.path.join(
	os.path.dirname(__file__), 'reference_data.csv')
#: File to write simulation data to
OUTPUT_DATA_PATH = os.path.join(
	os.path.dirname(__file__), 'test_output.csv')
#: Length of simulation to run for both generation and checking
DURATION = 50
#: Number of cells with which to initialize simulation
NUM_CELLS = 2
#: Seconds in a cell cycle before forcibly triggering division
SECONDS_TO_DIVISION = 15


def process_path_timeseries_for_csv(path_ts):
	# type: (dict) -> dict
	'''Prepare path timeseries data for writing to CSV

	Processing Steps:

	1. Convert any dictionary keys that are tuples to strings by joining
	   the tuple elements with ``-``.
	2. Remove from the timeseries any data where each timepoint is not
	   numeric. For example, we remove data where each timepoint is a
	   matrix or a string.

	.. note:: We assume that across timepoints for the same variable,
		the data types are the same.

	Returns:
		dict: A timeseries that can be saved to a CSV with
		:py:func:`save_flat_timeseries`.
	'''
	# Convert tuple keys to strings
	str_keys = dict()
	for key, value in path_ts.items():
		try:
			iter(key)
			if not isinstance(key, str):
				key = ",".join(key)
		except TypeError:
			pass
		str_keys[key] = value

	remove_keys = [
		# Remove matrices
		key for key, val in str_keys.items()
		if isinstance(val[0], list)
	]
	# Remove keys with non-numeric data
	for key, val in str_keys.items():
		try:
			float(val[0])
		except (ValueError, TypeError):
			remove_keys.append(key)
	for key in set(remove_keys):
		del str_keys[key]
	return str_keys


def save_flat_timeseries(
	timeseries,
	out_dir=os.path.dirname(REFERENCE_DATA_PATH),
	filename=os.path.basename(REFERENCE_DATA_PATH),
):
	# type: (dict, str, str) -> None
	'''Save a flattened timeseries to a CSV file

	The CSV file will have one column for each key in the timeseries.
	The heading will be the key, and the rows will contain the
	timeseries data, one row per timepoint, in increasing order of time.
	'''
	n_rows = max([len(val) for val in timeseries.values()])
	rows = [{} for _ in range(n_rows)]
	for key, val in timeseries.items():
		for i, elem in enumerate(val):
			rows[i][key] = elem
	with io.open(os.path.join(out_dir, filename), 'wb') as f:
		writer = tsv.dict_writer(f, timeseries.keys(), delimiter=',')
		writer.writeheader()
		writer.writerows(rows)


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
		emitter_config, DURATION, NUM_CELLS, SECONDS_TO_DIVISION)
	path_ts = emitter.get_path_timeseries()
	processed_ts = process_path_timeseries_for_csv(path_ts)
	return processed_ts


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'--generate', '-g',
		action='store_true',
		default=False,
		help='Save reference data to {}.'.format(REFERENCE_DATA_PATH),
	)
	parser.add_argument(
		'--check', '-c',
		action='store_true',
		default=False,
		help='Check that model behavior matches {}.'.format(REFERENCE_DATA_PATH),
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
		reference_ts = load_timeseries(REFERENCE_DATA_PATH)
		# We load simulation data from the file so that all
		# transformations by CSV reader/writer are applied
		simulation_ts = load_timeseries(OUTPUT_DATA_PATH)
		assert_timeseries_close(simulation_ts, reference_ts)
	if args.generate:
		save_flat_timeseries(processed_ts)


if __name__ == '__main__':
	main()
