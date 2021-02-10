#! /usr/bin/env python

"""
Script to run all options and models for NCA.  Need to specify options and
models of interest with the OPTIONS and MODELS variables.  Saves log files
to log/ and uses the default output location for NCA results.
"""

import argparse
import csv
import os
import re
import subprocess
from typing import List

from wholecell.utils import parallelization


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
OUTPUT_DIR = os.path.join(FILE_LOCATION, 'out')
os.makedirs(OUTPUT_DIR, exist_ok=True)
LOG_DIR = os.path.join(FILE_LOCATION, 'log')
os.makedirs(LOG_DIR, exist_ok=True)
FOLD_CHANGE_SCRIPT = os.path.join(FILE_LOCATION, 'fold_changes.py')
SUMMARY_FILE = os.path.join(OUTPUT_DIR, 'summary.tsv')

OPTIONS = [
	'average-seq',
	'linear',
	'global-expression',
	'split',
	'sigma-factors',
	]
MODELS = [
	'robust_nca',
	'constrained_nca',
	]


def solve_nca(model: str, options: List[str]):
	label = '-'.join([model] + [option.split('-')[0] for option in options])
	flags = ' '.join([f'--{option}' for option in options])
	cmd = f'{FOLD_CHANGE_SCRIPT} -f -l {label} -m {model} {flags}'

	print(f'Running: {cmd}')
	with open(os.path.join(LOG_DIR, f'{label}.log'), 'w') as f:
		subprocess.run(cmd.split(), stdout=f)

def compile_output():
	"""
	Compiles output log files to a spreadsheet.
	"""

	with open(SUMMARY_FILE, 'w') as f:
		fieldnames = [
			'Options',
			'Overall match', 'Overall total', 'Overall percent',
			'Negative match', 'Negative total', 'Negative percent',
			'Positive match', 'Positive total', 'Positive percent',
			'Error',
			'r, all', 'n, all',
			'r, consistent', 'n, consistent'
			]
		writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
		writer.writeheader()

		for log_file in sorted(os.listdir(LOG_DIR)):
			if not log_file.endswith('.log'):
				continue

			with open(os.path.join(LOG_DIR, log_file)) as log:
				data = {'Options': log_file}
				for line in log.readlines():
					if line.startswith('Overall matches:'):
						match, total = re.findall('([0-9]*)/([0-9]*)', line)[0]
						percent = re.findall('([0-9]*\.[0-9]*)%', line)[0]
						data['Overall match'] = int(match)
						data['Overall total'] = int(total)
						data['Overall percent'] = float(percent)
					elif line.startswith('Negative regulation matches:'):
						match, total = re.findall('([0-9]*)/([0-9]*)', line)[0]
						percent = re.findall('([0-9]*\.[0-9]*)%', line)[0]
						data['Negative match'] = int(match)
						data['Negative total'] = int(total)
						data['Negative percent'] = float(percent)
					elif line.startswith('Positive regulation matches:'):
						match, total = re.findall('([0-9]*)/([0-9]*)', line)[0]
						percent = re.findall('([0-9]*\.[0-9]*)%', line)[0]
						data['Positive match'] = int(match)
						data['Positive total'] = int(total)
						data['Positive percent'] = float(percent)
					elif line.startswith('Error fitting original data (excluding unpredicted samples):'):
						error = re.findall('([0-9]*\.[0-9]*)', line)[0]
						data['Error'] = float(error)
					elif line.startswith('All: r='):
						r = re.findall('r=([0-9]*\.[0-9]*)', line)[0]
						n = re.findall('n=([0-9]*)', line)[0]
						data['r, all'] = float(r)
						data['n, all'] = float(n)
					elif line.startswith('Only if WCM consistent: r='):
						r = re.findall('r=([0-9]*\.[0-9]*)', line)[0]
						n = re.findall('n=([0-9]*)', line)[0]
						data['r, consistent'] = float(r)
						data['n, consistent'] = float(n)

			writer.writerow(data)

def parse_args() -> argparse.Namespace:
	"""Parse command line args for options to run."""

	parser = argparse.ArgumentParser()

	# General options
	parser.add_argument('--compile',
		action='store_true',
		help='Compile log output into a spreadsheet without rerunning solvers (useful if compile_output() changed).')
	parser.add_argument('--cpus',
		type=int,
		default=2,
		help='Number of CPUs to use to run different option combinations. Beware of memory limitations with constrained_nca and many CPUs. (default: 2)')

	return parser.parse_args()

if __name__ == '__main__':
	args = parse_args()

	if args.compile:
		compile_output()
	else:
		all_args = [
			(model, [option for j, option in enumerate(OPTIONS) if i//2**j % 2 == 0])
			for model in MODELS
			for i in range(2**len(OPTIONS))
			]

		pool = parallelization.pool(args.cpus)
		results = [pool.apply_async(solve_nca, args) for args in all_args]
		pool.close()
		pool.join()

		compile_output()
