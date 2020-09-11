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
		raw_data = KnowledgeBaseEcoli()

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
	rnaIdToSymbol = {x['rna_id']: x['symbol'] for x in raw_data.genes}
	monomerToRna = {x['id'][:-3]: x['rna_id'][:-3] for x in sim_data.process.translation.monomer_data}
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
		output (str): path to tsv file with list of implemented genes
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
		output (str): path to tsv file with list of implemented genes
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
