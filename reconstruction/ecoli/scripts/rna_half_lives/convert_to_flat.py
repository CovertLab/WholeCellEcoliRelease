#! /usr/bin/env python
"""
Extracts half-lives of RNAs from Supplementary Table S4 of Chen et al.,
"Genome-wide study of mRNA degradation and transcript elongation in Escherichia
coli" (2014) and Supporting Table 5 of Bernstein et al., "Global analysis of
mRNA decay and abundance in Escherichi coli at single gene resoluation using
two-color fluorescent DNA microarrays" (2002). Optionally generates a comparison
plot for degradation rates reported for the same gene from these two sources.
"""

import ast
import io
import itertools
import os
import time
from typing import Dict, Set, Tuple

import numpy as np
import matplotlib.pyplot as plt

from wholecell.io import tsv
from wholecell.utils.filepath import ROOT_PATH

# If True, plot comparison between half-lives from Chen et al. and Bernstein
# et al.
PLOT_COMPARISON = False

# Number of decimals when rounding half-life values (in units of min)
ROUND_N_DECIMALS = 1

# Directories
FILE_LOCATION = os.path.realpath(os.path.dirname(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
OUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUT_DIR):
	os.makedirs(OUT_DIR)

# Files
CHEN_INPUT = os.path.join(DATA_DIR, 'msb145794-sup-0009-supp_table_s4.csv')  # Supplementary Table S4
BERNSTEIN_INPUT = os.path.join(DATA_DIR, '3181Table5.csv')  # Supporting Table 5, saved as .csv (originally .xls) and comments removed
GENES_FILE = os.path.join(
	ROOT_PATH, 'reconstruction', 'ecoli', 'flat', 'genes.tsv')
TRANSCRIPTION_UNITS_FILE = os.path.join(
	ROOT_PATH, 'reconstruction', 'ecoli', 'flat', 'transcription_units.tsv')
OUTPUT_FLAT_FILE = os.path.join(
	ROOT_PATH, 'reconstruction', 'ecoli', 'flat', 'rna_half_lives.tsv')
if PLOT_COMPARISON:
	OUTPUT_COMPARISON_PLOT_FILE = os.path.join(
		OUT_DIR, 'chen_bernstein_comparison.png')


def get_symbols_to_ids():
	# type: () -> Tuple[Dict[str, str], Set[str]]
	"""
	Builds a mapping from gene/transcription unit symbols to their corresponding
	IDs that gets used in later data import functions.

	Returns:
		symbols_to_ids: Dictionary that maps gene/transcription unit symbols
		to gene/transcription unit IDs used by the model.
		all_gene_ids: Set of all gene IDs used by the model.
	"""
	symbols_to_ids = {}
	ids_to_synonyms = {}
	all_gene_ids = set()
	all_gene_symbols = set()

	# Map gene symbols to gene IDs
	with io.open(GENES_FILE, 'rb') as f:
		reader = tsv.reader(f, delimiter='\t')

		headers = next(reader)
		while headers[0].startswith('#'):
			headers = next(reader)

		gene_id_index = headers.index('id')
		gene_symbol_index = headers.index('symbol')
		synonyms_index = headers.index('synonyms')

		for line in reader:
			gene_symbol = line[gene_symbol_index]
			gene_id = line[gene_id_index]
			synonyms = ast.literal_eval(line[synonyms_index])

			symbols_to_ids[gene_symbol] = gene_id
			ids_to_synonyms[gene_id] = synonyms
			all_gene_ids.add(gene_id)
			all_gene_symbols.add(gene_symbol)

	# Add nonconflicting gene ID synonyms to mapping (this is needed because
	# Bernstein et al. uses older genomic annotations)
	for (gene_id, synonyms) in ids_to_synonyms.items():
		for synonym in synonyms:
			if synonym not in symbols_to_ids:
				symbols_to_ids[synonym] = gene_id
				all_gene_symbols.add(gene_symbol)

	# Map transcription unit symbols to transcription unit IDs
	with io.open(TRANSCRIPTION_UNITS_FILE, 'rb') as f:
		reader = tsv.reader(f, delimiter='\t')

		headers = next(reader)
		while headers[0].startswith('#'):
			headers = next(reader)

		tu_id_index = headers.index('id')
		tu_symbol_index = headers.index('common_name')
		gene_list_index = headers.index('genes')

		for line in reader:
			tu_symbol = line[tu_symbol_index]
			gene_list = ast.literal_eval(line[gene_list_index])

			# Skip TUs whose symbols are just existing single gene symbols
			if tu_symbol not in all_gene_symbols:
				if len(gene_list) > 1:
					symbols_to_ids[line[tu_symbol_index]] = line[tu_id_index]
				else:
					symbols_to_ids[line[tu_symbol_index]] = gene_list[0]

	return symbols_to_ids, all_gene_ids


def get_chen_half_lives(symbols_to_ids):
	# type: (Dict[str, str]) -> Dict[str, float]
	"""
	Reads RNA half lives reported by Chen et al., and stores them in a
	dictionary. Any gene/transcription unit listed that does not exist in the
	model are skipped.

	Args:
		symbols_to_ids: Dictionary that maps gene/transcription unit symbols
		to gene/transcription unit IDs used by the model.

	Returns:
		id_to_half_life: Dictionary that maps gene/transcription unit IDs to
		their half lives (in minutes) reported by Chen et al.
	"""
	id_to_half_life = {}

	with io.open(CHEN_INPUT, 'rb') as f:
		reader = tsv.reader(f, delimiter=',')

		headers = next(reader)
		rna_symbol_index = headers.index('RNA')
		avg_lifetime_index = headers.index('Avg lifetime')

		for line in reader:
			rna_symbol = line[rna_symbol_index]
			avg_lifetime = float(line[avg_lifetime_index])

			if rna_symbol in symbols_to_ids:
				# Convert average liftime to half-life
				id_to_half_life[symbols_to_ids[rna_symbol]] = avg_lifetime * np.log(2)

	return id_to_half_life


def get_bernstein_half_lives(symbols_to_ids):
	# type: (Dict[str, str]) -> Dict[str, float]
	"""
	Reads RNA half lives reported by Bernstein et al., and stores them in a
	dictionary. Any gene listed that does not exist in the model are skipped.
	Bernstein et al. only reports half-lives of monocistronic RNAs.

	Args:
		symbols_to_ids: Dictionary that maps gene/transcription unit symbols
		to gene/transcription unit IDs used by the model.

	Returns:
		id_to_half_life: Dictionary that maps gene IDs to their RNA half lives
		(in minutes) reported by Bernstein et al.
	"""
	id_to_half_life = {}

	with io.open(BERNSTEIN_INPUT, 'rb') as f:
		reader = tsv.reader(f, delimiter=',')

		headers = next(reader)
		rna_symbol_index = headers.index('Symbol')
		half_life_index = headers.index('Half-lives in M9')

		for line in reader:
			if len(line[half_life_index]) == 0:
				continue

			rna_symbol = line[rna_symbol_index]
			half_life = float(line[half_life_index])

			if rna_symbol in symbols_to_ids:
				id_to_half_life[symbols_to_ids[rna_symbol]] = half_life

	return id_to_half_life


def build_half_life_table(chen_half_lives, bernstein_half_lives, all_gene_ids):
	# type: (Dict[str, float], Dict[str, float], Set[str]) -> None
	"""
	Builds the RNA half-life flat file that is used as raw_data for the
	simulation using the data from the two papers. If data for a gene exists in
	both papers, we use data from Chen et al. assuming that their methodology is
	more accurate. In the output file monocistronic RNAs (represented by their
	corresponding gene IDs) are listed first, followed by polycistronic RNAs.
	"""
	monocistronic_data = []
	polycistronic_data = []
	all_chen_gene_ids = set()

	for (rna_id, half_life) in chen_half_lives.items():
		if rna_id in all_gene_ids:
			monocistronic_data.append((rna_id, round(half_life, ROUND_N_DECIMALS)))
			all_chen_gene_ids.add(rna_id)
		else:
			polycistronic_data.append((rna_id, round(half_life, ROUND_N_DECIMALS)))

	# Only use data from Bernstein et al. if gene does not exist in Chen et al.
	for (gene_id, half_life) in bernstein_half_lives.items():
		if gene_id not in all_chen_gene_ids:
			monocistronic_data.append((gene_id, round(half_life, ROUND_N_DECIMALS)))

	# Sort by ID
	monocistronic_data.sort(key=lambda v: v[0])
	polycistronic_data.sort(key=lambda v: v[0])

	# Write to flat file
	with io.open(OUTPUT_FLAT_FILE, 'wb') as f:
		print('Writing to {}'.format(f.name))
		writer = tsv.writer(f, quotechar="'", lineterminator='\n')
		writer.writerow(['# Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(['"id"', '"half_life (units.min)"'])

		for row in itertools.chain(monocistronic_data, polycistronic_data):
			writer.writerow([f'"{row[0]}"', f'{row[1]}'])


def compare_half_lives(chen_half_lives, bernstein_half_lives):
	# type: (Dict[str, float], Dict[str, float]) -> None
	"""
	Generates a scatter plot that compares the half-lives of RNAs reported in
	both Chen et al. and Bernstein et al.
	"""
	chen_genes = set(chen_half_lives.keys())
	bernstein_genes = set(bernstein_half_lives.keys())

	overlapping_genes = sorted(list(chen_genes & bernstein_genes))

	overlapping_chen_half_lives = np.array(
		[chen_half_lives[gene] for gene in overlapping_genes])
	overlapping_bernstein_half_lives = np.array(
		[bernstein_half_lives[gene] for gene in overlapping_genes])

	plt.figure(figsize=(8, 8))
	plt.scatter(
		overlapping_bernstein_half_lives, overlapping_chen_half_lives, s=2)
	plt.xlim([0, 35])
	plt.ylim([0, 35])
	plt.xlabel('RNA half-lives from Bernstein et al. (min)')
	plt.ylabel('RNA half-lives from Chen et al. (min)')
	plt.savefig(OUTPUT_COMPARISON_PLOT_FILE)


if __name__ == '__main__':
	symbols_to_ids, all_gene_ids = get_symbols_to_ids()
	chen_half_lives = get_chen_half_lives(symbols_to_ids)
	bernstein_half_lives = get_bernstein_half_lives(symbols_to_ids)

	build_half_life_table(chen_half_lives, bernstein_half_lives, all_gene_ids)
	if PLOT_COMPARISON:
		compare_half_lives(chen_half_lives, bernstein_half_lives)
