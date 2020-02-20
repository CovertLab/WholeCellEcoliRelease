#! /usr/bin/env python

"""
Extracts average concentrations for metabolites from Lempp et al. Systematic
identification of metabolites controlling gene expression in E. coli. 2019.

Data in lempp2019.tsv was converted to a tsv from the file downloaded from
https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12474-1/MediaObjects/41467_2019_12474_MOESM4_ESM.xlsx
"""

from __future__ import absolute_import, division, print_function

import csv
import os
import urllib

import numpy as np
from typing import Dict, Iterable


FILE_LOCATION = os.path.realpath(os.path.dirname(__file__))
CONCENTRATIONS_FILE = os.path.join(FILE_LOCATION, 'lempp2019.tsv')
OUTPUT_FILE = os.path.join(FILE_LOCATION, 'lempp_concentrations.tsv')


def metabolite_concentrations():
	# type: () -> Dict[str, float]
	"""
	Load average metabolite concentrations at the first time point.

	Returns:
		met_conc: KEGG molecule ID to concentration (in M)
	"""

	met_conc = {}

	with open(CONCENTRATIONS_FILE) as f:
		reader = csv.reader(f, delimiter='\t')

		start_conc_col = reader.next().index('intracellular concentrations (\xc2\xb5M)')
		reader.next()  # discard line
		n_conc = np.sum([t.startswith('t0') for t in reader.next()[start_conc_col:]])
		end_conc_col = start_conc_col + n_conc
		id_col = reader.next().index('KEGG')

		for line in reader:
			met_id = line[id_col]
			try:
				conc = np.array(line[start_conc_col:end_conc_col], float).mean()
			except ValueError as e:
				# Concentration data does not exist ('-')
				continue

			# Convert from uM to M concentration
			met_conc[met_id] = conc / 1e6

	return met_conc

def kegg_to_ecocyc(kegg_ids):
	# type: (Iterable[str]) -> Dict[str, str]
	"""
	Retrieve EcoCyc ID for each KEGG ID.

	Args:
		kegg_ids: KEGG molecule IDs to get EcoCyc IDs for

	Returns:
		mapping: KEGG ID to EcoCyc ID
	"""

	mapping = {}
	id_type = 'Kegg:'
	url = 'https://websvc.biocyc.org/ECOLI/foreignid?ids='
	ids = ','.join(['{}{}'.format(id_type, i) for i in kegg_ids])

	u = urllib.urlopen(url + ids)
	reader = csv.reader(u, delimiter='\t')

	for line in reader:
		if line[1] == '1':
			mol_id = line[0].split(id_type)[1]
			ecocyc_id = line[2]
			mapping[mol_id] = ecocyc_id

	return mapping

def save_concentrations(conc, mapping):
	# type: (Dict[str, float], Dict[str, str]) -> None
	"""
	Save EcoCyc ID concentrations to a file.

	Args:
		conc: ID to concentration (in M)
		mapping: ID to EcoCyc ID
	"""

	conc = {mapping[m]: c for m, c in conc.items() if m in mapping}

	with open(OUTPUT_FILE, 'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Metabolite', 'Lempp Concentration (units.mol/units.L)'])

		for m, c in sorted(conc.items(), key=lambda d: d[0]):
			writer.writerow([m, '{:.2e}'.format(c).replace('e+0', 'e').replace('e-0', 'e-')])


if __name__ == '__main__':
	concentrations = metabolite_concentrations()
	kegg_mapping = kegg_to_ecocyc(concentrations.keys())
	save_concentrations(concentrations, kegg_mapping)
