#! /usr/bin/env python

"""
Inspects raw_data and sim_data to identify genes, metabolites and
kinetic constraints that are implemented in the model.

Will need to be updated with new processes and when submodels are extended.

Outputs:
	tsv files with lists of implemented features in the model
		functional_genes.tsv: genes that have a functional role
		metabolite_pools.tsv: metabolites that have a concentration
		kinetic_constraints.tsv: metabolic reaction kinetic constraints
		regulation.tsv: transcriptional regulation implemented
"""

from __future__ import absolute_import, division, print_function

import argparse
import io
import os
import time

import numpy as np
from six.moves import cPickle, zip

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
from wholecell.io import tsv
from wholecell.utils import filepath


FILE_LOCATION = os.path.dirname(os.path.realpath(__file__))
GENES_FILE = 'functional_genes.tsv'
METABOLITES_FILE = 'metabolite_pools.tsv'
KINETICS_FILE = 'kinetic_constraints.tsv'
REGULATION_FILE = 'regulation.tsv'
AMINO_ACIDS_FILE_BASE = 'amino_acid'


def load_raw_data(path):
	"""
	Loads raw_data from an existing file, otherwise loads from flat files.

	Args:
		path (str): path to raw_data cPickle file to load
	"""

	if os.path.exists(path):
		print('Loading raw_data from {}...'.format(path))
		with open(path, 'rb') as f:
			raw_data = cPickle.load(f)
	else:
		print('Loading raw_data from flat files...')
		raw_data = KnowledgeBaseEcoli(operons_on=False)

	return raw_data

def load_sim_data(path, raw_data):
	"""
	Loads sim_data from an existing file, otherwise recalculates it.

	Args:
		path (str): path to sim_data cPickle file to load
		raw_data (KnowledgeBaseEcoli object): raw data for simulations
	"""

	if os.path.exists(path):
		print('Loading sim_data from {}...'.format(path))
		with open(path, 'rb') as f:
			sim_data = cPickle.load(f)
	else:
		print('Calculating sim_data...')
		sim_data = fitSimData_1(raw_data)

	return sim_data

def save_genes(raw_data, sim_data, output):
	"""
	Gets genes that are implemented in the model and saves the list to a tsv file.

	Args:
		raw_data (KnowledgeBaseEcoli object): raw data for simulations
		sim_data (SimulationDataEcoli object): simulation data
		output (str): path to tsv file with list of implemented genes
	"""

	validMonomers = sim_data.process.translation.monomer_data['id']

	monomers = []

	##### Metabolism #####
	metMonomers = [
		x for x in sim_data.process.metabolism.catalyst_ids
		if x not in sim_data.process.complexation.complex_names
		and x not in sim_data.process.equilibrium.complex_name_to_rxn_idx
		]
	metComplexes = [
		x for x in sim_data.process.metabolism.catalyst_ids
		if x in sim_data.process.complexation.complex_names
		or x in sim_data.process.equilibrium.complex_name_to_rxn_idx
		]

	assert len(metMonomers) + len(metComplexes) == len(sim_data.process.metabolism.catalyst_ids)

	for metComplex in metComplexes:
		if metComplex in sim_data.process.complexation.complex_names:
			metMonomers += sim_data.process.complexation.get_monomers(metComplex)['subunitIds'].tolist()
		elif metComplex in sim_data.process.equilibrium.complex_name_to_rxn_idx:
			for subunit in sim_data.process.equilibrium.get_monomers(metComplex)['subunitIds'].tolist():
				if subunit in sim_data.process.complexation.complex_names:
					metMonomers += sim_data.process.complexation.get_monomers(subunit)['subunitIds'].tolist()
				elif subunit in validMonomers:
					metMonomers += [subunit]
		else:
			raise Exception

	monomers += metMonomers

	##### Translation #####
	translationMonomers = []
	translationMonomers += sim_data.process.complexation.get_monomers(sim_data.molecule_ids.s30_full_complex)['subunitIds'].tolist()
	translationMonomers += sim_data.process.complexation.get_monomers(sim_data.molecule_ids.s50_full_complex)['subunitIds'].tolist()

	monomers += translationMonomers

	##### Transcription #####
	transcriptionMonomers = []
	transcriptionMonomers += sim_data.process.complexation.get_monomers(sim_data.molecule_ids.full_RNAP)['subunitIds'].tolist()

	monomers += transcriptionMonomers

	##### RNA Decay #####
	rnaDecayMonomers = []
	rnaDecayMonomers += sim_data.process.rna_decay.endoRNase_ids
	rnaDecayMonomers += sim_data.molecule_groups.exoRNases

	monomers += rnaDecayMonomers

	##### Transcriptional Regulation #####
	tfMonomers = []
	tfIds = [x + '[c]' for x in sim_data.process.transcription_regulation.tf_to_tf_type]
	tfComplexes = [
		x for x in tfIds
		if x in sim_data.process.complexation.complex_names
		or x in sim_data.process.equilibrium.complex_name_to_rxn_idx
		or x in sim_data.process.two_component_system.complex_to_monomer
		]
	tfMonomers += [
		x for x in tfIds
		if x not in sim_data.process.complexation.complex_names
		and x not in sim_data.process.equilibrium.complex_name_to_rxn_idx
		and x not in sim_data.process.two_component_system.complex_to_monomer
		]

	assert len(tfMonomers) + len(tfComplexes) == len(tfIds)

	for tfComplex in tfComplexes:
		if tfComplex in sim_data.process.complexation.complex_names:
			tfMonomers += sim_data.process.complexation.get_monomers(tfComplex)['subunitIds'].tolist()
		elif tfComplex in sim_data.process.equilibrium.complex_name_to_rxn_idx:
			for subunit in sim_data.process.equilibrium.get_monomers(tfComplex)['subunitIds'].tolist():
				if subunit in sim_data.process.complexation.complex_names:
					tfMonomers += sim_data.process.complexation.get_monomers(subunit)['subunitIds'].tolist()
				elif subunit in validMonomers:
					tfMonomers += [subunit]
		elif tfComplex in sim_data.process.two_component_system.complex_to_monomer:
			for subunit in sim_data.process.two_component_system.complex_to_monomer[tfComplex]:
				if subunit in sim_data.process.complexation.complex_names:
					tfMonomers += sim_data.process.complexation.get_monomers(subunit)['subunitIds'].tolist()
				elif subunit in validMonomers:
					tfMonomers += [subunit]
		else:
			raise Exception

	monomers += tfMonomers

	monomers = set([x for x in monomers if x in validMonomers])

	# Get gene names for each monomer implemented
	rnaIdToSymbol = {x['rna_ids'][0]: x['symbol'] for x in raw_data.genes}
	monomerToRna = {x['id'][:-3]: x['cistron_id'] for x in sim_data.process.translation.monomer_data}
	geneNames = [rnaIdToSymbol[monomerToRna[monomer[:-3]]] for monomer in monomers]

	# Save data to output tsv file
	## Update with new processes
	functional_monomers = {
		'Metabolism': metMonomers,
		'Translation': translationMonomers,
		'Transcription': transcriptionMonomers,
		'RNA Decay': rnaDecayMonomers,
		'Transcription Regulation': tfMonomers,
		}

	with io.open(output, 'wb') as f:
		print('\nWriting gene info to {}'.format(output))
		writer = tsv.writer(f)
		writer.writerow(['Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(['Gene', 'RNA ID', 'Monomer ID', 'Process'])

		monomers_added = set()
		for gene, monomer in sorted(zip(geneNames, monomers), key=lambda v: v[0]):
			processes = []
			for function, subset in functional_monomers.items():
				if monomer in subset:
					processes.append(function)

			# Prevent duplicates in different compartments
			monomer = monomer[:-3]
			if monomer not in monomers_added:
				writer.writerow([gene, monomerToRna[monomer], monomer, ', '.join(processes)])
				monomers_added.add(monomer)

	print('Number of genes: {}'.format(len(monomers_added)))

def save_metabolites(raw_data, sim_data, output):
	"""
	Gets metabolites that have a concentration in the model and saves the list to
	a tsv file.

	Args:
		raw_data (KnowledgeBaseEcoli object): raw data for simulations
		sim_data (SimulationDataEcoli object): simulation data
		output (str): path to tsv file with list of metabolites with a concentration
	"""

	metabolites = list(sim_data.process.metabolism.conc_dict.keys())
	bennett = [m['Metabolite'] for m in raw_data.metabolite_concentrations
		if not np.isnan(m['Bennett Concentration'].asNumber())]
	lempp = [m['Metabolite'] for m in raw_data.metabolite_concentrations
		if not np.isnan(m['Lempp Concentration'].asNumber())]

	with io.open(output, 'wb') as f:
		print('\nWriting metabolite info to {}'.format(output))
		writer = tsv.writer(f)
		writer.writerow(['Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(['Metabolite', 'Source'])

		for m in sorted(metabolites):
			sources = []
			if m[:-3] in bennett:
				sources.append('Bennett et al. 2009')
			if m[:-3] in lempp:
				sources.append('Lempp et al. 2019')

			if sources:
				source = ', '.join(sources)
			else:
				source = 'Biomass (EcoCyc GEM, Bremer and Dennis. 1996., and Neidhardt. 2006.)'
			writer.writerow([m, source])

	print('Number of metabolites: {}'.format(len(metabolites)))

def save_kinetics(raw_data, sim_data, output):
	"""
	Gets kinetic constraints that are used in the model and saves the list to
	a tsv file.

	Args:
		raw_data (KnowledgeBaseEcoli object): raw data
		sim_data (SimulationDataEcoli object): simulation data
		output (str): path to tsv file with list of kinetic constraints
	"""

	kinetic_constraints = sim_data.process.metabolism.extract_kinetic_constraints(
		raw_data, sim_data)

	with io.open(output, 'wb') as f:
		print('\nWriting kinetics info to {}'.format(output))
		writer = tsv.writer(f)
		writer.writerow(['Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(['Reaction ID', 'Enzyme', 'Adjusted kcat (1/s)', 'Saturation'])

		# TODO: include info like pubmed ID, substrates, KM, KI, temp
		for (rxn, enz), c in sorted(kinetic_constraints.items(), key=lambda d: d[0][0]):
			writer.writerow([rxn, enz, c['kcat'], c['saturation']])

	print('Number of kinetic constraints: {}'.format(len(kinetic_constraints)))

def save_amino_acid_pathways(sim_data, output):
	"""
	Saves parameters related to mechanistic amino acid biosynthesis and transport
	to a tsv file.

	Args:
		sim_data (SimulationDataEcoli object): simulation data
		output (str): path to a file to save
	"""

	metabolism = sim_data.process.metabolism
	aa_ids = [aa[:-3] for aa in sim_data.molecule_groups.amino_acids]

	# Synthesis
	synthesis_output = output + '_synthesis.tsv'
	with open(synthesis_output, 'wb') as f:
		print('\nWriting amino acid synthesis info to {}'.format(synthesis_output))
		writer = tsv.writer(f)
		writer.writerow(['Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(['Amino acid', 'kcat forward (1/s)',
			'kcat reverse or degradation (1/s)', 'KI (uM)', 'Upstream KM (mM)',
			'Reverse KM (mM)', 'Degradation KM (mM)',
			])
		for aa, kcat_fwd, kcat_rev, ki, upstream_km, reverse_km, deg_km in zip(
				aa_ids,
				metabolism.aa_kcats_fwd,
				metabolism.aa_kcats_rev,
				metabolism.aa_kis,
				metabolism.aa_upstream_kms,
				metabolism.aa_reverse_kms,
				metabolism.aa_degradation_kms,
				):
			# Not functinoally implemented and most parameters do not exist
			if aa == 'L-SELENOCYSTEINE':
				continue
			writer.writerow([aa, f'{kcat_fwd:.3g}', f'{kcat_rev:.3g}',
				f'{1e6 * ki:.3g}' if np.isfinite(ki) else '',
				', '.join([f'{a}: {k:.3g}' for a, k in zip(aa_ids, 1e3 * upstream_km) if k > 0]),
				f'{1e3 * reverse_km:.3g}' if np.isfinite(reverse_km) else '',
				f'{1e3 * deg_km:.3g}' if np.isfinite(deg_km) else '',
				])

	# Transport
	transport_output = output + '_transport.tsv'
	with open(transport_output, 'wb') as f:
		print('\nWriting amino acid transport info to {}'.format(transport_output))
		writer = tsv.writer(f)
		writer.writerow(['Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(['Amino acid', 'kcat import (1/s)',
			'kcat export (1/s)', 'Import KI (mM)', 'Export KM (mM)',
			])
		for aa, kcat_im, kcat_ex, ki_im, km_ex in zip(
				aa_ids,
				metabolism.import_kcats_per_aa,
				metabolism.export_kcats_per_aa,
				metabolism.aa_import_kis,
				metabolism.aa_export_kms,
				):
			# Not functinoally implemented and most parameters do not exist
			if aa == 'L-SELENOCYSTEINE':
				continue
			writer.writerow([aa, f'{kcat_im:.3g}', f'{kcat_ex:.3g}',
				f'{1e3 * ki_im:.3g}' if np.isfinite(ki_im) else '',
				f'{1e3 * km_ex:.3g}' if np.isfinite(km_ex) else '',
				])

def save_regulation(sim_data, output):
	"""
	Gets transcription regulation incorporated in the model and saves the list to
	a tsv file.

	Args:
		sim_data (SimulationDataEcoli object): simulation data
		output (str): path to tsv file with list of regulation
	"""

	def count_regulation(indices, direction, regulators, cistrons):
		genes = [cistrons[i] for i in np.unique(indices)]
		total_pos_regulation = np.sum(direction == 1)
		total_neg_regulation = np.sum(direction == -1)

		cistron_regulation = {}
		for i, d in zip(indices, direction):
			cistron_regulation.setdefault(i, set()).add(d)
		n_genes_both = 0
		n_genes_pos = 0
		n_genes_neg = 0
		regulation_type = []
		for i, d in sorted(cistron_regulation.items(), key=lambda d: d[0]):
			if len(d) == 2:
				regulation_type.append('both')
				n_genes_both += 1
			elif 1 in d:
				regulation_type.append('+')
				n_genes_pos += 1
			else:
				regulation_type.append('-')
				n_genes_neg += 1

		regulation = dict(zip(genes, regulation_type))
		n_regulators = len(np.unique(regulators))

		return regulation, total_pos_regulation, total_neg_regulation, n_genes_both, n_genes_pos, n_genes_neg, n_regulators

	transcription = sim_data.process.transcription
	replication = sim_data.process.replication
	delta_prob = sim_data.process.transcription_regulation.delta_prob

	cistron_ids = transcription.cistron_data['id']
	cistron_id_to_idx = {
		cistron: i
		for i, cistron in enumerate(cistron_ids)
		}
	cistron_id_to_symbol = {
		gene['cistron_id']: gene['symbol']
		for gene in replication.gene_data
		}
	attenuation_k = transcription.attenuation_k.asNumber()

	# Cistron indices for regulation
	tf_regulated_indices = np.array([
		i for i, v in zip(delta_prob['deltaI'], delta_prob['deltaV'])
		if v != 0
		])
	ppgpp_regulated_indices = np.array([
		cistron_id_to_idx[cistron]
		for cistron in transcription.ppgpp_regulated_genes
		])
	atten_regulated_indices = transcription.attenuated_rna_indices[np.where(attenuation_k)[1]]

	# Direction for regulation
	tf_regulated_direction = np.sign(delta_prob['deltaV'][delta_prob['deltaV'] != 0])
	ppgpp_regulated_direction = np.sign(transcription.ppgpp_fold_changes)
	atten_regulated_direction = np.sign(attenuation_k[attenuation_k != 0])

	# Regulator indices
	tf_regulator_indices = np.array([
		j for j, v in zip(delta_prob['deltaJ'], delta_prob['deltaV'])
		if v != 0
		])
	ppgpp_regulator_indices = np.zeros_like(ppgpp_regulated_indices)  # Only regulated by ppGpp
	atten_regulator_indices = np.where(attenuation_k)[0]

	# Combined data
	all_regulated_indices = np.hstack((
		tf_regulated_indices,
		ppgpp_regulated_indices,
		atten_regulated_indices,
	))
	all_regulated_direction = np.hstack((
		tf_regulated_direction,
		ppgpp_regulated_direction,
		atten_regulated_direction,
	))
	all_regulator_indices = np.hstack((
		[f'tf-{i}' for i in tf_regulator_indices],
		[f'ppgpp-{i}' for i in ppgpp_regulator_indices],
		[f'atten-{i}' for i in atten_regulator_indices],
	))

	data = {}
	data['All regulation'] = count_regulation(all_regulated_indices,
		all_regulated_direction, all_regulator_indices, cistron_ids)
	data['TF regulation'] = count_regulation(tf_regulated_indices,
		tf_regulated_direction, tf_regulator_indices, cistron_ids)
	data['ppGpp regulation'] = count_regulation(ppgpp_regulated_indices,
		ppgpp_regulated_direction, ppgpp_regulator_indices, cistron_ids)
	data['Transcription attenuation'] = count_regulation(atten_regulated_indices,
		atten_regulated_direction, atten_regulator_indices, cistron_ids)

	all_regulated_cistrons = sorted(data['All regulation'][0])

	with io.open(output, 'wb') as f:
		print('\nWriting regulation info to {}'.format(output))
		writer = tsv.writer(f)
		writer.writerow(['Generated by {} on {}'.format(__file__, time.ctime())])
		for label, d in data.items():
			writer.writerow([
				f'{label} has {d[1] + d[2]} total regulatory interactions'
				f' ({d[1]} positive, {d[2]} negative)'
				f' with {d[3] + d[4] + d[5]} genes regulated'
				f' ({d[3]} positive and negative, {d[4]} positive, {d[5]} negative)'
				f' from {d[6]} regulators'
				])
		writer.writerow(['Cistron', 'Symbol'] + list(data.keys()))

		for cistron in all_regulated_cistrons:
			writer.writerow([cistron, cistron_id_to_symbol[cistron]] + [d[0].get(cistron, 'None') for d in data.values()])

	print(f'Number of genes regulated: {len(data["All regulation"][0])}')
	print(f'Number of regulation interactions: {data["All regulation"][1] + data["All regulation"][2]}')

def parse_args():
	"""
	Parses command line arguments.

	Returns:
		argparse.Namespace object: parsed arguments and values
	"""

	default_output_dir = FILE_LOCATION

	parser = argparse.ArgumentParser(description='Script to save lists of'
		' included genes, metabolites and kinetic constraints in the model')

	parser.add_argument('-r', '--raw-data', default='',
		help='Path to raw_data cPickle object to load, recalculates raw_data if not specified')
	parser.add_argument('-s', '--sim-data', default='',
		help='Path to sim_data cPickle object to load, recalculates sim_data if not specified')
	parser.add_argument('-o', '--output', default=default_output_dir,
		help='Directory path to save tsv files (default: {})'.format(default_output_dir))

	return parser.parse_args()


if __name__ == '__main__':
	# Parse command line args
	args = parse_args()

	# Load required data
	raw_data = load_raw_data(args.raw_data)
	sim_data = load_sim_data(args.sim_data, raw_data)

	# Analyze data and save tsv files
	filepath.makedirs(args.output)
	save_genes(raw_data, sim_data, os.path.join(args.output, GENES_FILE))
	save_metabolites(raw_data, sim_data, os.path.join(args.output, METABOLITES_FILE))
	save_kinetics(raw_data, sim_data, os.path.join(args.output, KINETICS_FILE))
	save_amino_acid_pathways(sim_data, os.path.join(args.output, AMINO_ACIDS_FILE_BASE))
	save_regulation(sim_data, os.path.join(args.output, REGULATION_FILE))
