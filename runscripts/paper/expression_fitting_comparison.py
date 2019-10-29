"""

Extracts information about the differences between gene expression, with and
without fitting parameter values (that is, for ribosome and RNA polymerase
expression).

Note that these are hard paths to simulation output - these will need to be
modified to reevaluate this script on new data.

"""

from __future__ import absolute_import, division, print_function

import cPickle
import os
import csv

import numpy as np

# ID parsing constants

RRNA_ID_TO_GENE = {
	'RRSA-RRNA[c]': 'rrsA',
	'RRFA-RRNA[c]': 'rrfA',
	'RRLA-RRNA[c]': 'rrlA'
	}

TRNA_ID_TO_GENE = {
	'RNA0-300[c]': 'valZ',
	'RNA0-301[c]': 'lysY',
	'RNA0-302[c]': 'lysZ',
	'RNA0-303[c]': 'lysQ',
	'RNA0-304[c]': 'asnW',
	'RNA0-305[c]': 'ileY',
	'RNA0-306[c]': 'metV'
	}

SLICE_RNA_TO_FRAME_ID = slice(None, -len('_RNA[c]'))
TRNA_TRAILING_ID = '-tRNA[c]'
SLICE_TRNA_TO_GENE_NAME = slice(None, -len(TRNA_TRAILING_ID))

# Output configuration

MAX_DECIMALS = 3 # precision of changes we care about
MAX_GENE_PRINT = 15 # more than this, and we just print a number rather than names

# Printing

PRINTED_DIVIDER_SIZE = 79


def main(unfit_path, fit_path, condition = 'basal'):
	# Load sim data (parameters)

	unfit_parameters = cPickle.load(open(UNFIT_SOURCE))
	fit_parameters = cPickle.load(open(FIT_SOURCE))

	# Extract data and parse names

	# NOTE: This should really assert that these parameters (names, etc.) are the
	# same between the two files

	rna_data = unfit_parameters.process.transcription.rnaData.fullArray()

	rna_ids = rna_data['id']

	get_subunits = lambda names: unfit_parameters.process.complexation.getMonomers(
		names
		)['subunitIds']

	ribosome30S_ids = get_subunits(unfit_parameters.moleculeIds.s30_fullComplex)
	ribosome50S_ids = get_subunits(unfit_parameters.moleculeIds.s50_fullComplex)

	rnap_ids = get_subunits(unfit_parameters.moleculeIds.rnapFull)

	protein_to_rna = {
		entry['id']: entry['rnaId']
		for entry in unfit_parameters.process.translation.monomerData.fullArray()
		}

	rna_to_index = {
		entry['id']: i
		for (i, entry)
		in enumerate(rna_data)
		}

	r_rna_ids = []
	r_protein_ids = []

	for ids in (ribosome30S_ids, ribosome50S_ids):
		for id_ in ids:
			try:
				r_protein_ids.append(protein_to_rna[id_])

			except KeyError:
				r_rna_ids.append(id_)

	rnap_ids = [protein_to_rna[id_] for id_ in rnap_ids]

	trna_ids = rna_ids[rna_data['isTRna']]

	id_to_name = {}

	# Parse IDs to common gene names

	# Load the frame-to-gene mapping

	# This file should eventually be available elsewhere, but had to be extracted
	# from another branch and commit.
	# (053211cf57ffcfe77436607345a1963e91ee94b5/models/ecoli/analysis/causal_network/Genes.txt)
	with open(os.path.join('runscripts', 'paper', 'Genes.txt')) as f:
		for row in csv.DictReader(f, dialect = 'excel-tab'):
			id_to_name[row['Gene Name']] = row['Names'].split(' // ')[0].strip('"')

	# Pick out the different groupings
	r_rna_indices = [rna_to_index[id_] for id_ in r_rna_ids]
	r_protein_indices = [rna_to_index[id_] for id_ in r_protein_ids]
	rnap_indices = [rna_to_index[id_] for id_ in rnap_ids]
	trna_indices = [rna_to_index[id_] for id_ in trna_ids]

	is_else = np.ones(rna_data.size, np.bool)

	is_else[r_rna_indices] = False
	is_else[r_protein_indices] = False
	is_else[rnap_indices] = False
	is_else[trna_indices] = False

	# Collect gene names (parsing from IDs)
	r_rna_names = [
		RRNA_ID_TO_GENE[id_]
		for id_ in r_rna_ids
		]
	r_protein_names = [
		id_to_name[id_[SLICE_RNA_TO_FRAME_ID]]
		for id_ in r_protein_ids
		]
	rnap_names = [
		id_to_name[id_[SLICE_RNA_TO_FRAME_ID]]
		for id_ in rnap_ids
		]
	trna_names = [
		(
			id_[SLICE_TRNA_TO_GENE_NAME]
			if id_.endswith(TRNA_TRAILING_ID)
			else TRNA_ID_TO_GENE[id_]
			)
		for id_ in trna_ids
		]
	else_names = [
		(
			id_to_name[id_[SLICE_RNA_TO_FRAME_ID]]
			if (id_[SLICE_RNA_TO_FRAME_ID] in id_to_name)
			else id_[SLICE_RNA_TO_FRAME_ID]
			)
		for id_ in rna_ids[is_else]
		]

	# Calculations

	# Get synthesis probabilities
	unfit_synth = unfit_parameters.process.transcription.rnaSynthProb[condition]
	fit_synth = fit_parameters.process.transcription.rnaSynthProb[condition]

	# Compute ratios
	ratios = fit_synth/unfit_synth

	ratios = np.round(ratios, MAX_DECIMALS) # Clip to desired significance

	# Make sure everything that was zero (no expression), stayed zero, and
	# vice-versa (since otherwise ratios are hard to interpret)
	became_nonzero = (fit_synth != 0) & (unfit_synth == 0)
	became_zero = (fit_synth == 0) & (unfit_synth != 0)

	assert not became_nonzero.any()
	assert not became_zero.any()

	# Output tables of fold-changes
	def print_table(entry_names, ratios, cutoff = MAX_GENE_PRINT):
		entry_names = np.array(entry_names)

		isfinite = np.where(np.isfinite(ratios))

		entry_names = entry_names[isfinite]
		ratios = ratios[isfinite]

		sorting = np.argsort(entry_names)

		entry_names = entry_names[sorting]
		ratios = ratios[sorting]

		(values, inverse, counts) = np.unique(
			ratios,
			return_inverse = True,
			return_counts = True
			)

		elements = [np.where(inverse == i)[0] for i in xrange(values.size)]

		for (r, e, c) in zip(values, elements, counts):
			if c < cutoff:
				names = ', '.join(entry_names[e])

			else:
				names = '{} genes'.format(c)

			print('{names}\t{value:0.{precision}F}'.format(
				names = names,
				value = r,
				precision = MAX_DECIMALS
				))

	print('rProteins')
	print('-'*PRINTED_DIVIDER_SIZE)
	print_table(r_protein_names, ratios[r_protein_indices])

	print()
	print('RNA polymerase subunits')
	print('-'*PRINTED_DIVIDER_SIZE)
	print_table(rnap_names, ratios[rnap_indices])

	print()
	print('rRNAs')
	print('-'*PRINTED_DIVIDER_SIZE)
	print_table(r_rna_names, ratios[r_rna_indices])

	print()
	print('tRNAs')
	print('-'*PRINTED_DIVIDER_SIZE)
	print_table(trna_names, ratios[trna_indices])

	print()
	print('All other proteins')
	print('-'*PRINTED_DIVIDER_SIZE)
	print_table(else_names, ratios[is_else])

if __name__ == '__main__':
	# Paths
	SET_ROOT = os.path.join('/', 'scratch', 'PI', 'mcovert', 'wc_ecoli', 'paper')
	PATH_TO_SIMDATA = os.path.join('kb', 'simData_Most_Fit.cPickle')
	UNFIT_SOURCE = os.path.join(
		SET_ROOT,
		'SET_L',
		'20190910.155143__SET_L_4_gens_256_seeds_3_conditions_unfit_ribosome_and_rna_poly_expression',
		PATH_TO_SIMDATA
	)
	FIT_SOURCE = os.path.join(
		SET_ROOT,
		'SET_C',
		'20190906.151026__SET_C_4_gens_256_seeds_3_conditions_with_growth_noise_and_D_period',
		PATH_TO_SIMDATA
	)

	for condition in ['basal', 'with_aa', 'no_oxygen']:
		print()
		print('='*PRINTED_DIVIDER_SIZE)
		print('Condition: {}'.format(condition))
		print()
		main(UNFIT_SOURCE, FIT_SOURCE, condition)
