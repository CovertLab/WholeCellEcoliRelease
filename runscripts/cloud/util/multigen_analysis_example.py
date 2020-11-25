"""A small example for CustomWorkflow to run within a Docker Firetask.
Given a seed dir, run example multigen analysis on its sim generations.
"""

import os
import sys

import numpy as np

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader


def print_time_range(sim_out_dir):
	# type: (str) -> float
	"""Print the time range from the given cell simOut."""
	times = TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time")
	print(f'    Time range: [{times[0]} .. {times[-1]}]')
	return times[-1]


def print_bulk_count_range(sim_out_dir):
	# type: (str) -> np.ndarray
	"""Print the start and end BulkMolecules counts from the given cell simOut."""
	counts = TableReader(os.path.join(sim_out_dir, "BulkMolecules")).readColumn("counts")
	print(f'    BulkMolecules counts start / end: {counts[0]} / {counts[-1]}')
	return counts[-1]


def analyze(seed_dir):
	# type: (str) -> None
	"""Examine each simOut dir within the seed_dir."""
	paths = AnalysisPaths(seed_dir, multi_gen_plot=True)
	last_time = 0.0
	last_counts = np.arange(0)

	# Print some data from each simOut dir.
	for cell_dir in paths.get_cells():
		print(cell_dir)
		sim_out_dir = os.path.join(cell_dir, 'simOut', '')
		last_time = print_time_range(sim_out_dir)
		last_counts = print_bulk_count_range(sim_out_dir)

	# Write an output file into seed_dir/../count_out/<seed_key>/complex/
	base_dir, seed_key = os.path.split(os.path.dirname(os.path.join(seed_dir, '')))
	output_dir = os.path.join(base_dir, 'count_out', seed_key, 'complex')
	os.makedirs(output_dir, exist_ok=True)
	output_file = os.path.join(output_dir, 'multigen_analysis_example.txt')

	with open(output_file, 'w') as out:
		print(f'{last_time=}', file=out)
		print(f'{last_counts=}', file=out)
		print(f'Wrote {output_file}')


if __name__ == '__main__':
	seed_dir_arg = sys.argv[1]  # e.g. /wcEcoli/out/counts/wildtype_000000/000000
	analyze(seed_dir_arg)
