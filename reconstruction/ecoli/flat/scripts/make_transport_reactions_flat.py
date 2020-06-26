'''
Makes a tsv file with all transport reaction ids

goes through all reactions in reaction.tsv and compares substrate and product ids

'''

from __future__ import absolute_import, division, print_function

import csv
import io
import os
import re
from typing import Dict

import six

from reconstruction.spreadsheets import read_tsv
from wholecell.io import tsv


TSV_DIALECT = csv.excel_tab

REACTIONS_FILE = os.path.join("reconstruction", "ecoli", "flat", "reactions.tsv")
OUT_FILE = os.path.join("out", "transport_reactions.tsv")

# make list of transport reactions
transport_reactions = []
for row in read_tsv(REACTIONS_FILE):
	reaction_id = row["reaction id"]  # type: str
	stoichiometry = row["stoichiometry"]  # type: Dict[str, int]

	# get substrates and products
	substrates = [mol_id for mol_id, coeff in six.viewitems(stoichiometry) if coeff < 0]
	products = [mol_id for mol_id, coeff in six.viewitems(stoichiometry) if coeff > 0]
	substrates_no_loc = [re.sub("[[@*&?].*[]@*&?]", "", mol_id) for mol_id in substrates]
	products_no_loc = [re.sub("[[@*&?].*[]@*&?]", "", mol_id) for mol_id in products]

	overlap_no_loc = set(substrates_no_loc) & set(products_no_loc)

	# if overlap between substrate and product names with no location:
	for mol_id in list(overlap_no_loc):
		sub = [mol for mol in substrates if mol_id in mol]
		prod = [mol for mol in products if mol_id in mol]
		overlap = set(sub) & set(prod)

		# if there is no overlap between those substrates and products with locations included
		if len(overlap) == 0:
			# print('sub ' + str(sub))
			# print('prod ' + str(prod))
			transport_reactions.append(
				reaction_id.encode('ascii', 'ignore').decode('ascii'))

# sort reactions to save them in ordered list
transport_reactions = list(set(transport_reactions))
transport_reactions = sorted(transport_reactions)

# save list of transport reactions
if os.path.exists(OUT_FILE):
	os.remove(OUT_FILE)

with io.open(OUT_FILE, 'ab') as tsvfile:
	writer = tsv.writer(tsvfile, quoting=csv.QUOTE_NONNUMERIC)
	writer.writerow(["reaction id"])
	for reaction_id_ in transport_reactions:
		append_line = [reaction_id_]
		writer.writerow(append_line)
