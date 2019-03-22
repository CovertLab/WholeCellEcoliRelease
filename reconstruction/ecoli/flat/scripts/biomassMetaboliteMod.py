'''
Modify metabolites file so that it no longer stores biomass data in a stupid way
'''
import numpy as np
from os.path import isfile, join
import csv
import collections
from reconstruction.spreadsheets import JsonReader, JsonWriter
CSV_DIALECT = csv.excel_tab
FLAT_DIR = '/home/users/nruggero/Repos/wcEcoli/reconstruction/ecoli/flat'

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

# Rebuild metabolites biomass columns
metaboliteData = load_tsv(join(FLAT_DIR,'metabolites.tsv'))
for metabolite in metaboliteData:
	core_concentration = []
	core_location = []
	for data in  metabolite['biomassInfo']['core']:
		core_concentration.append(data['mmol/gDCW'])
		core_location.append(data['location'])
	wildtype_concentration = []
	wildtype_location = []

	for data in  metabolite['biomassInfo']['wildtype']:
		wildtype_concentration.append(data['mmol/gDCW'])
		wildtype_location.append(data['location'])

	metabolite.pop('biomassInfo')
	metabolite['core_concentration'] = core_concentration
	metabolite['core_location'] = core_location
	metabolite['wildtype_concentration'] = wildtype_concentration
	metabolite['wildtype_location'] = wildtype_location

write_tsv(join(FLAT_DIR,'metabolites.tsv'), metaboliteData)
