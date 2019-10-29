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

import argparse
import cPickle
import csv
import os
import time

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.ecoli.fit_sim_data_1 import fitSimData_1
from wholecell.utils import filepath, units


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

	validMonomers = sim_data.process.translation.monomerData['id']

	monomers = []

	##### Metabolism #####
	metMonomers = [
		x for x in sim_data.process.metabolism.catalystsList
		if x not in sim_data.process.complexation.complexNames
		and x not in sim_data.process.equilibrium.complexNameToRxnIdx
		]
	metComplexes = [
		x for x in sim_data.process.metabolism.catalystsList
		if x in sim_data.process.complexation.complexNames
		or x in sim_data.process.equilibrium.complexNameToRxnIdx
		]

	assert len(metMonomers) + len(metComplexes) == len(sim_data.process.metabolism.catalystsList)

	for metComplex in metComplexes:
		if metComplex in sim_data.process.complexation.complexNames:
			metMonomers += sim_data.process.complexation.getMonomers(metComplex)['subunitIds'].tolist()
		elif metComplex in sim_data.process.equilibrium.complexNameToRxnIdx:
			for subunit in sim_data.process.equilibrium.getMonomers(metComplex)['subunitIds'].tolist():
				if subunit in sim_data.process.complexation.complexNames:
					metMonomers += sim_data.process.complexation.getMonomers(subunit)['subunitIds'].tolist()
				elif subunit in validMonomers:
					metMonomers += [subunit]
		else:
			raise Exception

	monomers += metMonomers

	##### Translation #####
	translationMonomers = []
	translationMonomers += sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)['subunitIds'].tolist()
	translationMonomers += sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)['subunitIds'].tolist()

	monomers += translationMonomers

	##### Transcription #####
	transcriptionMonomers = []
	transcriptionMonomers += sim_data.process.complexation.getMonomers(sim_data.moleculeIds.rnapFull)['subunitIds'].tolist()

	monomers += transcriptionMonomers

	##### RNA Decay #####
	rnaDecayMonomers = []
	rnaDecayMonomers += sim_data.process.rna_decay.endoRnaseIds
	rnaDecayMonomers += sim_data.moleculeGroups.exoRnaseIds

	monomers += rnaDecayMonomers

	##### Transcriptional Regulation #####
	tfMonomers = []
	tfIds = [x + '[c]' for x in sim_data.process.transcription_regulation.tfToTfType]
	tfComplexes = [
		x for x in tfIds
		if x in sim_data.process.complexation.complexNames
		or x in sim_data.process.equilibrium.complexNameToRxnIdx
		or x in sim_data.process.two_component_system.complexToMonomer
		]
	tfMonomers += [
		x for x in tfIds
		if x not in sim_data.process.complexation.complexNames
		and x not in sim_data.process.equilibrium.complexNameToRxnIdx
		and x not in sim_data.process.two_component_system.complexToMonomer
		]

	assert len(tfMonomers) + len(tfComplexes) == len(tfIds)

	for tfComplex in tfComplexes:
		if tfComplex in sim_data.process.complexation.complexNames:
			tfMonomers += sim_data.process.complexation.getMonomers(tfComplex)['subunitIds'].tolist()
		elif tfComplex in sim_data.process.equilibrium.complexNameToRxnIdx:
			for subunit in sim_data.process.equilibrium.getMonomers(tfComplex)['subunitIds'].tolist():
				if subunit in sim_data.process.complexation.complexNames:
					tfMonomers += sim_data.process.complexation.getMonomers(subunit)['subunitIds'].tolist()
				elif subunit in validMonomers:
					tfMonomers += [subunit]
		elif tfComplex in sim_data.process.two_component_system.complexToMonomer:
			for subunit in sim_data.process.two_component_system.complexToMonomer[tfComplex]:
				if subunit in sim_data.process.complexation.complexNames:
					tfMonomers += sim_data.process.complexation.getMonomers(subunit)['subunitIds'].tolist()
				elif subunit in validMonomers:
					tfMonomers += [subunit]
		else:
			raise Exception

	monomers += tfMonomers

	monomers = set([x for x in monomers if x in validMonomers])

	# Get gene names for each monomer implemented
	rnaIdToSymbol = {x['rnaId']: x['symbol'] for x in raw_data.genes}
	monomerToRna = {x['id'][:-3]: x['rnaId'][:-3] for x in sim_data.process.translation.monomerData}
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

	with open(output, 'w') as f:
		print('\nWriting gene info to {}'.format(output))
		writer = csv.writer(f, delimiter='\t')
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

	print 'Number of genes: {}'.format(len(monomers_added))

def save_metabolites(raw_data, sim_data, output):
	"""
	Gets metabolites that have a concentration in the model and saves the list to
	a tsv file.

	Args:
		raw_data (KnowledgeBaseEcoli object): raw data for simulations
		sim_data (SimulationDataEcoli object): simulation data
		output (str): path to tsv file with list of implemented genes
	"""

	metabolites = sim_data.process.metabolism.concDict.keys()
	bennett = [m['Metabolite'] for m in raw_data.metaboliteConcentrations]

	with open(output, 'w') as f:
		print('\nWriting metabolite info to {}'.format(output))
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(['Metabolite', 'Source'])

		for m in sorted(metabolites):
			if m[:-3] in bennett:
				source = 'Bennett et al. 2009'
			else:
				source = 'Biomass (EcoCyc GEM, Bremer and Dennis. 1996., and Neidhardt. 2006.)'
			writer.writerow([m, source])

	print 'Number of metabolites: {}'.format(len(metabolites))

def save_kinetics(sim_data, output):
	"""
	Gets kinetic constraints that are used in the model and saves the list to
	a tsv file.

	Args:
		sim_data (SimulationDataEcoli object): simulation data
		output (str): path to tsv file with list of implemented genes
	"""

	kinetic_constraints = sim_data.process.metabolism.constraintDict
	disabled_constraints = set(sim_data.process.metabolism.constraintsToDisable)

	with open(output, 'w') as f:
		print('\nWriting kinetics info to {}'.format(output))
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Generated by {} on {}'.format(__file__, time.ctime())])
		writer.writerow(['Pubmed ID', 'Reaction ID', 'kcat (1/s)',
			'KM (uM)', 'Adjusted kcat (1/s)', 'Enzyme',
			'Substrates', 'Substrates for KM', 'Temperature (C)', 'Excluded'])

		for _, c in sorted(kinetic_constraints.items(), key=lambda d: d[0]):
			rxn = c['reactionID']
			excluded = rxn in disabled_constraints
			writer.writerow([c['Pubmed ID'], rxn,
				c['kcat'].asNumber(1 / units.s),
				c['kM'], '{:.4g}'.format(c['kcatAdjusted'].asNumber(1 / units.s)),
				c['enzymeIDs'], c['substrateIDs'], c['Concentration Substrates'],
				c['Temp'], excluded])

	print 'Number of kinetic constraints: {}'.format(len(kinetic_constraints))

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
	save_kinetics(sim_data, os.path.join(args.output, KINETICS_FILE))
