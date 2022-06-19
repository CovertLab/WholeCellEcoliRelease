#! /usr/bin/env python

"""
Convert WCM output tables to tsv files for easy portability.  Modify the DATA
variable to select table and columns you would like to convert.  You must supply
the path to the simOut directory with tables as the first and only arg to this
script.

TODO (Travis):
	process all tables/columns?
	integrate with script base class for easier selection of variant, seed, gen etc.
"""

import csv
import os
import sys

import numpy as np

from wholecell.io.tablereader import TableReader, SUBCOLUMNS_KEY


DATA = [
	('MonomerCounts', 'monomerCounts'),
	('FBAResults', 'reactionFluxes'),
	('EnzymeKinetics', 'metaboliteCountsFinal'),
	('EnzymeKinetics', 'targetFluxes'),
	('RnaSynthProb', 'rnaSynthProb'),
	('RnapData', 'rnaInitEvent')
	]


def get_labels(reader, column):
	subcolumns = reader.readAttribute(SUBCOLUMNS_KEY)
	return np.array(reader.readAttribute(subcolumns[column]))

def save_data(out_dir, t, data, labels, column):
	filename = os.path.join(out_dir, f'{column}.tsv')
	print(f'Converting to {filename}')
	with open(filename, 'w') as f:
		writer = csv.writer(f, delimiter='\t')

		writer.writerow(['Time (s)'] + list(t))
		for label, d in zip(labels, data.T):
			writer.writerow([label] + list(d))


if __name__ == '__main__':
	out_dir = sys.argv[1]

	main = TableReader(os.path.join(out_dir, 'Main'))
	t = main.readColumn('time')

	for table, column in DATA:
		reader = TableReader(os.path.join(out_dir, table))
		data = reader.readColumn(column)
		labels = get_labels(reader, column)

		save_data(out_dir, t, data, labels, column)
