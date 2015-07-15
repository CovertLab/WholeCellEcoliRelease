'''
Unify files that are used to create bulkMolecules state
'''
import numpy as np
from os.path import isfile, join
import csv
import collections
from reconstruction.spreadsheets import JsonReader, JsonWriter
CSV_DIALECT = csv.excel_tab
FLAT_DIR = 'reconstruction/ecoli/flat/'

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
	outfile = open(file_name, "w")
	writer = JsonWriter(outfile, fieldnames, dialect = CSV_DIALECT)
	writer.writeheader()
	writer.writerows(l)

chromosomeIds = ['CHROM_FORWARD', 'CHROM_REVERSE', 'CHROM_FORWARD_COMPLEMENT', 'CHROM_REVERSE_COMPLEMENT']
chromosomeLocations = ['c','c','c','c']
chromosomeMWs = [309648605.409, 307349976.4, 308237613.853, 308854963.746]

chromosomeData = []
for idx in range(len(chromosomeIds)):
	chrom = {
			"id" : chromosomeIds[idx],
			"location" : [chromosomeLocations[idx]],
			"mw" : [0.] * len(molecular_weight_order),
			}
	chrom['mw'][molecular_weight_order['DNA']] = chromosomeMWs[idx]
	chromosomeData.append(chrom)

write_tsv(join(FLAT_DIR,'chromosome.tsv'), chromosomeData)