#! /usr/bin/env python

"""
Extracts average concentrations for metabolites from Lempp et al. Systematic
identification of metabolites controlling gene expression in E. coli. 2019.

Data in lempp2019.tsv was converted to a tsv from the file downloaded from
https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-12474-1/MediaObjects/41467_2019_12474_MOESM4_ESM.xlsx

Data in park2016.tsv was converted to a tsv from supplementary table 5
(manually removed some columns and rows) in the file downloaded from
https://static-content.springer.com/esm/art%3A10.1038%2Fnchembio.2077/MediaObjects/41589_2016_BFnchembio2077_MOESM583_ESM.pdf

Data in kochanowski2017*.tsv was converted to a tsv from EV table 6-1
(absolute) and EV table 6-2 (relative) (manually removed some rows)
from the file downloaded from
https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20167402&file=msb167402-sup-0007-TableEV6.xlsx
"""

from __future__ import absolute_import, division, print_function

import csv
import os
import sys
import time
import urllib

import numpy as np
from typing import Any, Dict


# Directories
FILE_LOCATION = os.path.realpath(os.path.dirname(__file__))
DATA_DIR = os.path.join(FILE_LOCATION, 'data')
OUT_DIR = os.path.join(FILE_LOCATION, 'out')
if not os.path.exists(OUT_DIR):
	os.makedirs(OUT_DIR)

# Files
LEMPP_INPUT = os.path.join(DATA_DIR, 'lempp2019.tsv')
PARK_INPUT = os.path.join(DATA_DIR, 'park2016.tsv')
KOCHANOWSKI_ABSOLUTE_INPUT = os.path.join(DATA_DIR, 'kochanowski2017absolute.tsv')
KOCHANOWSKI_RELATIVE_INPUT = os.path.join(DATA_DIR, 'kochanowski2017relative.tsv')
ABSOLUTE_OUTPUT_FILE = os.path.join(OUT_DIR, '{}_concentrations.tsv')
RELATIVE_OUTPUT_FILE = os.path.join(OUT_DIR, 'relative_metabolite_concentrations.tsv')

# Correct EcoCyc IDs to match the whole-cell model ID
ECOCYC_SUBSTITUTIONS = {
	'D-GLUCOSAMINE-6-P': 'CPD-13469',
	'D-glucopyranose-6-phosphate': 'GLC-6-P',
	'Isocitrate': 'THREO-DS-ISO-CITRATE',
	'CPD-18719': 'FRUCTOSE-6P',
	}

# Kochanowski mappings
## Special media IDs used later (implemented in WCM)
GLC_MEDIA = 'minimal'
ACETATE_MEDIA = 'minimal_acetate'
SUCCINATE_MEDIA = 'minimal_succinate'
## Map metabolite column to wcm IDs (not all could be matched)
KOCHANOWSKI_METABOLITES = {
	'Glycerol-P': 'GLYCEROL-3P',
	'F1P': 'FRU1P',
	# 'Ga6P': '',
	'G6P': 'GLC-6-P',
	'F6P': 'FRUCTOSE-6P',
	'FBP': 'FRUCTOSE-16-DIPHOSPHATE',
	'DHAP': 'DIHYDROXY-ACETONE-PHOSPHATE',
	'BPG': 'DPG',
	# 'xPG': '',
	'PEP': 'PHOSPHO-ENOL-PYRUVATE',
	'6PG': 'CPD-2961',
	'Ru5P': 'RIBULOSE-5P',
	'R5P': 'RIBOSE-5P',
	'Xu5P': 'XYLULOSE-5-PHOSPHATE',
	'R1P': 'RIBOSE-1P',
	'S7P': 'D-SEDOHEPTULOSE-7-P',
	'Lactate': 'Lactate',
	'Acetyl-CoA': 'ACETYL-COA',
	'Citrate/Isocitrate': 'THREO-DS-ISO-CITRATE',
	'Aconitate': 'CIS-ACONITATE',
	'Alpha ketoglutarate': '2-KETOGLUTARATE',
	'Succinate': 'SUC',
	'Fumarate': 'FUM',
	'Malate': 'MAL',
	'AMP': 'AMP',
	'ADP': 'ADP',
	'ATP': 'ATP',
	'GMP': 'GMP',
	'GDP': 'GDP',
	'GTP': 'GTP',
	'IMP': 'IMP',
	'NAD': 'NAD',
	'NADH': 'NADH',
	'NADP': 'NADP',
	'NADPH': 'NADPH',
	'GTTred': 'GLUTATHIONE',
	'GTTox': 'OXIDIZED-GLUTATHIONE',
	'Asparagine': 'ASN',
	'Aspartate': 'L-ASPARTATE',
	'Arginine': 'ARG',
	'Glutamine': 'GLN',
	'Glutamate': 'GLT',
	'Phenylalanine': 'PHE',
	'Tyrosine': 'TYR',
	'Panthothenate': 'PANTOTHENATE',
	# 'UDP-hexose': '',
	'cAMP': 'CAMP',
	}
## Map media headers to wcm IDs (not all are currently valid media)
KOCHANOWSKI_MEDIA = {
	' M9 galactose 2g per L': 'minimal_galactose',
	' M9 acetate 3.6g per L': ACETATE_MEDIA,
	' M9 mannose 2g per L': 'minimal_mannose',
	' M9 pyruvate 5g per L': 'minimal_pyruvate',
	' M9 lactate 5g per L': 'minimal_lactate',
	' M9 glycerol 2g per L': 'minimal_glycerol',
	' M9 sorbitol 2g per L': 'minimal_sorbitol',
	' M9 fructose 2g per L': 'minimal_fructose',
	' M9 succinate 2g per L': SUCCINATE_MEDIA,
	' M9 glcNAc 2g per L': 'minimal_glcNAc',
	' M9 mannitol 2g per L': 'minimal_mannitol',
	' M9 gluconate 2g per L': 'minimal_gluconate',
	' M9 glucose 2g per L first experiment': GLC_MEDIA,  # methods say 5 g/L
	' M9 G6P 2g per L': 'minimal_g6p',
	' M9 glucose 2g per L plus CAA 2g per L': 'minimal_plus_cas_amino_acids',
	# ' M9 glucose 0 mueM Chloramphenicol': '',
	# ' M9 glucose 0.75 mueM Chloramphenicol': '',
	# ' M9 glucose 1 mueM Chloramphenicol': '',
	# ' M9 glucose 2 mueM Chloramphenicol': '',
	# ' M9 glucose 4 mueM Chloramphenicol': '',
	# ' M9 glucose 6 mueM Chloramphenicol': '',
	# ' M9 glucose 8 mueM Chloramphenicol': '',
	# ' M9 glucose 10 mueM Chloramphenicol': '',
	}


def lempp_concentrations():
	# type: () -> Dict[str, float]
	"""
	Load Lempp data for average metabolite concentrations at the first time point.

	Returns:
		met_conc: EcoCyc molecule ID to concentration (in M)
	"""

	met_conc = {}

	with open(LEMPP_INPUT) as f:
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
			except ValueError as _:
				# Concentration data does not exist ('-')
				continue

			# Convert from uM to M concentration
			met_conc[met_id] = conc / 1e6

	return kegg_to_ecocyc(met_conc)

def park_concentrations():
	# type: () -> Dict[str, float]
	"""
	Load Park data for reported metabolite concentrations.

	Returns:
		met_conc: EcoCyc molecule ID to concentration (in M)
	"""

	met_conc = {}

	with open(PARK_INPUT) as f:
		reader = csv.reader(f, delimiter='\t')

		reader.next()  # discard line
		headers = reader.next()
		id_col = headers.index('KEGG ID')
		conc_col = headers.index('E. coli')

		for line in reader:
			met_id = line[id_col]
			try:
				conc = float(line[conc_col])
			except ValueError as e:
				# Concentration data does not exist ('-')
				continue

			met_conc[met_id] = conc

	return kegg_to_ecocyc(met_conc)

def load_kochanowski(filename):
	# type: (str) -> (Dict[str, np.ndarray[float]], np.ndarray[str])
	"""
	Load Kochanowski data (absolute or relative).

	Args:
		filename: path to file to load

	Returns:
		met_conc: WCM ID to concentration (absolute or relative)
	"""

	met_conc = {}

	with open(filename) as f:
		reader = csv.reader(f, delimiter='\t')

		reader.next()  # discard line
		headers = reader.next()[1:]
		valid_conditions = np.array([h in KOCHANOWSKI_MEDIA for h in headers])
		condition_headers = np.array([KOCHANOWSKI_MEDIA.get(h) for h in headers])[valid_conditions]
		met_col = 0

		for line in reader:
			met_id = KOCHANOWSKI_METABOLITES.get(line[met_col])
			if met_id is None:
				continue

			conc = np.array(line[1:], float)[valid_conditions]
			met_conc[met_id] = conc
			conc[conc <= 0] = np.nan

	return met_conc, condition_headers

def kochanowski_concentrations():
	# type: () -> Dict[str, float]
	"""
	Load absolute Kochanowski concentration data in the glucose condition.

	Returns:
		met_conc: WCM ID to concentration (in M)
	"""

	raw_data, headers = load_kochanowski(KOCHANOWSKI_ABSOLUTE_INPUT)
	glc_col = np.where(headers == GLC_MEDIA)[0][0]

	met_conc = {}
	for met, conc in raw_data.items():
		conc_in_glc = conc[glc_col]
		if np.isfinite(conc_in_glc):
			met_conc[met] = conc_in_glc / 1000  # convert mM to M

	return met_conc

def kegg_to_ecocyc(data):
	# type: (Dict[str, Any]) -> Dict[str, Any]
	"""
	Convert a dictionary with KEGG ID keys to a dictionary with EcoCyc ID keys.

	Args:
		data: dictionary to convert KEGG molecule IDs to EcoCyc IDs

	Returns:
		new_data: new dictionary with EcoCyc IDs
	"""

	kegg_ids = data.keys()
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
			mapping[mol_id] = ECOCYC_SUBSTITUTIONS.get(ecocyc_id, ecocyc_id)

	new_data = {mapping[m]: c for m, c in data.items() if m in mapping}

	return new_data

def save_concentrations(conc, label):
	# type: (Dict[str, float], str) -> None
	"""
	Save EcoCyc ID concentrations to a file.

	Args:
		conc: ID to concentration (in M)
		label: column and filename dataset label (author)
	"""

	output = ABSOLUTE_OUTPUT_FILE.format(label.lower())
	with open(output, 'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Metabolite', '{} Concentration (units.mol/units.L)'.format(label)])

		for m, c in sorted(conc.items(), key=lambda d: d[0]):
			writer.writerow([m, '{:.2e}'.format(c).replace('e+0', 'e').replace('e-0', 'e-')])

def save_kochanowski_relative_changes():
	"""
	Convert relative Kochanowski concentration data to wcm IDs for all conditions.
	"""

	met_conc, headers = load_kochanowski(KOCHANOWSKI_RELATIVE_INPUT)

	# Reorder output to put most interesting conditions first
	first_headers = [GLC_MEDIA, ACETATE_MEDIA, SUCCINATE_MEDIA]
	reordered_headers = first_headers + [h for h in headers if h not in first_headers]
	reordered_indexing = np.array([np.where(headers == h)[0][0] for h in reordered_headers])

	with open(RELATIVE_OUTPUT_FILE, 'w') as f:
		writer = csv.writer(f, delimiter='\t', quotechar="'", lineterminator='\n')
		writer.writerow(['# Created with {} on {}'.format(' '.join(sys.argv), time.ctime())])
		writer.writerow(['Metabolite'] + reordered_headers)

		for met, changes in sorted(met_conc.items()):
			data = [
				change if np.isfinite(change) else 'NaN'
				for change in changes[reordered_indexing]
			]
			writer.writerow(['"{}"'.format(met)] + data)


if __name__ == '__main__':
	# Lempp 2019
	lempp = lempp_concentrations()
	save_concentrations(lempp, 'Lempp')

	# Park 2016
	park = park_concentrations()
	save_concentrations(park, 'Park')

	# Kochanowski 2017
	kochanowski = kochanowski_concentrations()
	save_concentrations(kochanowski, 'Kochanowski')
	save_kochanowski_relative_changes()
