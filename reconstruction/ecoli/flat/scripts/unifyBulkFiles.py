'''
Unify files that are used to create bulkMolecules state
'''
import numpy as np
from os.path import isfile, join
import csv
import collections
from reconstruction.spreadsheets import JsonReader, JsonWriter
CSV_DIALECT = csv.excel_tab
FLAT_DIR = '/home/users/nruggero/Repos/wcEcoli/reconstruction/ecoli/flat'

molecular_weight_keys = [
	'23srRNA',
	'16srRNA',
	'5srRNA',
	'tRNA',
	'mRNA',
	'miscRNA',
	'protein',
	'metabolite',
	'water',
	'DNA',
	'RNA' # nonspecific RNA
	]

molecular_weight_order = collections.OrderedDict([
			(key, index) for index, key in enumerate(molecular_weight_keys)
			])

def load_tsv(file_name):
	data = []

	with open(file_name) as csvfile:
		reader = JsonReader(csvfile, dialect = CSV_DIALECT)
		for row in reader:
			data.append(dict([(x, y) for x,y in row.iteritems()]))
	return data

def write_tsv(file_name, l):
	fieldnames = l[0].keys()
	outfile = open(join(file_name, file_name), "w")
	writer = JsonWriter(outfile, fieldnames, dialect = CSV_DIALECT)
	writer.writeheader()
	writer.writerows(l)

# Load all compartments
compartmentData = load_tsv(join(FLAT_DIR,'compartments.tsv'))
allCompartmentAbbrev = [x['abbrev'] for x in compartmentData]

# Rebuild metabolites
metaboliteData = load_tsv(join(FLAT_DIR,'metabolites.tsv'))
for metabolite in metaboliteData:
	metabolite['location'] = allCompartmentAbbrev
	if metabolite['id'] != 'H2O':
		metabolite['mw'] = [0.] * len(molecular_weight_order)
		metabolite['mw'][molecular_weight_order['metabolite']] = metabolite['mw7.2']
	elif metabolite['id'] == 'H2O':
		metabolite['mw'] = [0.] * len(molecular_weight_order)
		metabolite['mw'][molecular_weight_order['water']] = metabolite['mw7.2']
	else:
		raise Exception, 'Messed up masses!'
	metabolite.pop('mediaConc')
	metabolite.pop('maxExchange')
write_tsv(join(FLAT_DIR,'metabolites.tsv'), metaboliteData)

# Rebuild proteins
proteinData = load_tsv(join(FLAT_DIR,'proteins.tsv'))
for protein in proteinData:
	protein['location'] = [protein['location']]
	mwSave = protein['mw']
	protein['mw'] = [0.] * len(molecular_weight_order)
	protein['mw'][molecular_weight_order['protein']] = mwSave
write_tsv(join(FLAT_DIR,'proteins.tsv'), proteinData)

# Rebuild RNAs
rnaData = load_tsv(join(FLAT_DIR,'rnas.tsv'))
for rna in rnaData:
	rna['location'] = [rna['location']]
write_tsv(join(FLAT_DIR,'rnas.tsv'), rnaData)

# Rebuild protein complexes
proteinComplexData = load_tsv(join(FLAT_DIR,'proteinComplexes.tsv'))
for proteinComplex in proteinComplexData:
	proteinComplex['location'] = [proteinComplex['location']]
write_tsv(join(FLAT_DIR,'proteinComplexes.tsv'), proteinComplexData)

# Rebuild polymerized
polymerizedData = load_tsv(join(FLAT_DIR,'polymerized.tsv'))
for polymerized in polymerizedData:
	polymerized['location'] = allCompartmentAbbrev
	mwSave = polymerized['mw']
	polymerized['mw'] = [0.] * len(molecular_weight_order)
	polymerized['mw'][polymerized['mass key']] = mwSave
	polymerized.pop('mass key')
write_tsv(join(FLAT_DIR,'polymerized.tsv'), polymerizedData)