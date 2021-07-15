#! /usr/bin/env python

"""
Aggregate data from Javi's repo to compare to fold_changes.tsv included in wcm.
"""

from __future__ import absolute_import, division, print_function

import argparse
import io
import os
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np
from six.moves import zip

from reconstruction.spreadsheets import tsv_writer
from wholecell.io import tsv


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
SRC_FILE = os.path.join(FILE_LOCATION, 'kb.tsv')
WCM_FILE = os.path.join(FILE_LOCATION, 'fold_changes.tsv')
SHIFTS_FILE = os.path.join(FILE_LOCATION, 'shifts.tsv')
GENES_FILES = os.path.join(FILE_LOCATION, 'gene_names.tsv')
OUTPUT_FILE = os.path.join(FILE_LOCATION, 'updated_fold_changes.tsv')


def load_file(filename):
	# type: (str) -> List[List[str]]
	"""Load a tsv file."""

	with io.open(filename, 'rb') as f:
		reader = tsv.reader(f)  # type: Iterable
		return list(reader)

def load_src(attempt_match):
	# type: (bool) -> Dict[str, Dict[str, Dict[str, float]]]
	"""
	Loads and extracts gene regulation data from the source data.

	Args:
		attempt_match: if True, handles data in way that was most likely the
			original processing, otherwise uses expected directionality

	Returns:
		data: mean and standard deviation for each regulatory pair
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}
	"""

	raw_data = load_file(SRC_FILE)
	shifts = load_shifts()

	# Extract fold changes from data
	data = {}  # type: Dict[str, Dict[str, List[float]]]
	for line in raw_data:
		# Columns of interest
		condition = (line[4], line[5])
		tf = line[7]
		regulated_gene = line[9]
		magnitude = float(line[10])

		if tf not in data:
			data[tf] = {}
		direction = shifts[condition][tf]

		data[tf][regulated_gene] = data[tf].get(regulated_gene, []) + [direction * magnitude]

	# Calculate mean and std from data
	processed_data = {}
	for tf, regulated in data.items():
		tf_data = {}
		for gene, fcs in regulated.items():
			if attempt_match:
				fcs1 = np.abs(fcs)  # Bad!! - ignores annotated condition comparison regulation direction
			else:
				fcs1 = np.array(fcs)
			tf_data[gene] = {
				'mean': float(np.mean(fcs1)),
				'std': float(np.std(fcs1, ddof=1)),
				}
		processed_data[tf] = tf_data

	return processed_data

def load_wcm(attempt_match):
	# type: (bool) -> Dict[str, Dict[str, Dict[str, float]]]
	"""
	Loads and extracts gene regulation data from the whole-cell model.

	Args:
		attempt_match: if True, handles data in way that was most likely the
			original processing, otherwise uses expected directionality

	Returns:
		data: mean and standard deviation for each regulatory pair
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}
	"""

	raw_data = load_file(WCM_FILE)[1:]

	# Extract mean and std from data
	data = {}  # type: Dict[str, Dict[str, Dict[str, Any]]]
	for line in raw_data:
		if line[0].startswith('#'):
			continue

		# Columns of interest
		tf = line[0].strip()
		regulated_gene = line[1].strip()
		mean = float(line[2])
		std = float(line[3])
		if attempt_match:
			sign = 1
		else:
			sign = int(np.sign(float(line[5])))

		uncertain = float(line[5]) > 2

		if tf not in data:
			data[tf] = {}

		data[tf][regulated_gene] = {
			'mean': sign * mean,
			'std': std,
			'uncertain': uncertain,
			}

	return data

def load_shifts():
	# type: () -> Dict[Tuple[str, str], Dict[str, int]]
	"""
	Load TF regulation information for each shift.

	Returns:
		shifts: gene regulation direction for each condition
			1 for active TF in condition1 vs condition2
			-1 for active TF in condition2 vs condition1
			{(condition1, condition2): {TF1: 1, TF2: -1, ...}}
	"""

	genes = {int(line[0]): line[2] for line in load_file(GENES_FILES)}
	shifts = {}
	for line in load_file(SHIFTS_FILE):
		condition = (line[1], line[2])

		gene_dir = {genes[np.abs(int(v))]: int(np.sign(int(v))) for v in line[12:26] if v != '0'}
		shifts[condition] = gene_dir

	return shifts

def compare_data(data1, data2, desc, verbose=True):
	# type: (Dict[str, Dict[str, Dict[str, float]]], Dict[str, Dict[str, Dict[str, float]]], str, bool) -> None
	"""
	Compares regulation from source data to wcm data.

	Args:
		data1: mean and standard deviation for each regulatory pair
		data2: mean and standard deviation for each regulatory pair
		desc: description for comparison
		verbose: if True, prints additional regulation information

	Notes:
		dictionary structure for both inputs:
			{TF gene name: {regulated gene: {'mean': mean, 'std': std}}}
	"""

	# Track statistics
	total_regulation = {}
	discrepancies = {}
	direction_discrepancies = {}
	no_match = {}

	print('\n*** {} ***'.format(desc))

	# Print regulation that is different in the datasets
	if verbose:
		print('TF -> regulated gene: data1 vs data2')
	for tf, regulation in data1.items():
		total_regulation[tf] = len(regulation)
		discrepancies[tf] = 0
		direction_discrepancies[tf] = 0
		no_match[tf] = 0

		for gene, data in regulation.items():
			mean1 = data['mean']
			mean2 = data2.get(tf, {}).get(gene, {}).get('mean', 0)
			if np.abs(mean1 - mean2) > 0.01 and mean2 != 0:
				if np.sign(mean1) != np.sign(mean2):
					direction_discrepancies[tf] += 1
				discrepancies[tf] += 1
				if verbose:
					uncertain = '\tuncertain' if data2.get(tf, {}).get(gene, {}).get('uncertain', False) else ''
					print('{} -> {}: {:.2f} vs {:.2f}{}'.format(tf, gene, mean1, mean2, uncertain))
			elif mean2 == 0:
				no_match[tf] += 1
	if verbose:
		print()

	# Print summary statistics for each TF
	total_direction_discrepancies = 0
	total_discrepancies = 0
	total_no_match = 0
	total_interactions = 0
	print('TF: opposite direction, different mean, not in second dataset, total')
	for tf, total in total_regulation.items():
		tf_dir_disc = direction_discrepancies[tf]
		tf_disc = discrepancies[tf]
		tf_no_match = no_match[tf]

		total_direction_discrepancies += tf_dir_disc
		total_discrepancies += tf_disc
		total_no_match += tf_no_match
		total_interactions += total

		if verbose:
			print('{:5s}: {:3} {:3} {:3} {:3}'.format(tf, tf_dir_disc, tf_disc, tf_no_match, total))
	print('Total: {:3} {:3} {:3} {:3}'.format(total_direction_discrepancies, total_discrepancies, total_no_match, total_interactions))

def replace_uncertain_entries(data):
	# type: (Dict[str, Dict[str, Dict[str, float]]]) -> None
	"""
	Replaces entries that were marked as uncertain in fold changes
	(Regulation_direct 3 or 4) if the mean fold change is the opposite sign
	when correcting for the regulation direction.

	Args:
		data: mean and standard deviation for each regulatory pair
	"""

	raw_data = load_file(WCM_FILE)
	headers = raw_data[0]

	with tsv_writer(OUTPUT_FILE, headers) as writer:
		# Check each fold change and update if needed
		for line in raw_data[1:]:  # type: List[str]
			tf = line[0].strip('#" ')
			gene = line[1].strip()
			direction = int(line[5])

			# Update value if uncertain fold change and opposite direction as mean
			if np.abs(direction) > 2:  # uncertain if category 3 or 4
				entry = data.get(tf, {}).get(gene, {})
				mean = entry.get('mean', 0)
				std = entry.get('std', 0)
				new_direction = int(np.sign(mean))
				if new_direction * direction < 0:  # both have opposite signs
					line[2] = '{:.2f}'.format(np.abs(mean))
					line[3] = '{:.2f}'.format(std)
					line[5] = '{:.0f}'.format(new_direction)

			# Write updated lines to file
			row = [tf, gene, float(line[2]), float(line[3]), float(line[4]),
				int(line[5]), int(line[6]), int(line[7]), float(line[8])]
			d = {header: value for header, value in zip(headers, row)}
			writer.writerow(d)

def parse_args():
	# type: () -> argparse.Namespace
	"""
	Parses arguments from the command line.

	Returns:
		values of variables parsed from the command line
	"""

	parser = argparse.ArgumentParser()

	parser.add_argument('-m', '--match',
		action='store_true',
		help='If set, processes source data to best match wcm.')
	parser.add_argument('-v', '--verbose',
		action='store_true',
		help='If set, prints more information.')
	parser.add_argument('--replace',
		action='store_true',
		help='Replace uncertain fold changes if new direction is different.')

	return parser.parse_args()

if __name__ == '__main__':
	args = parse_args()

	src_data = load_src(args.match)
	wcm_data = load_wcm(args.match)

	compare_data(src_data, wcm_data, 'Javi repo vs wcm', args.verbose)

	if args.replace:
		replace_uncertain_entries(src_data)
