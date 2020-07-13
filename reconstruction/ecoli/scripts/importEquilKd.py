#!/usr/bin/env python

# Generates a tsv output file to replace equilibriumReactions.tsv with updated Kd values
# Requires 2 files: validation/ecoli/flat/conformationTF.tsv and reconstruction/ecoli/flat/equilibriumReactions.tsv
# Outputs TFoutput.tsv to directory that script is run from
# Need to replace ' with " in file after script runs

from __future__ import absolute_import, division, print_function

import csv
import io
import os
from typing import Any, Dict

import numpy as np
import matplotlib
matplotlib.use("Agg")

from wholecell.io import tsv


DATA_FILE = os.path.join("validation","ecoli","flat","conformationTF.tsv")
REACTION_FILE = os.path.join("reconstruction","ecoli","flat","equilibriumReactions.tsv")
OUTPUT_FILE = "TFoutput.tsv"

reactionDict = {}
with io.open(DATA_FILE, "rb") as csvfile:
	reader = tsv.dict_reader(csvfile)
	for row in reader:  # type: Dict[str, Any]
		if row["<Kd> (uM)"] != '' and row["<Kd> (uM)"] != '?':
			for reaction in row["EcoCyc ID reaction (metabolite vs. TF)"].split(", "):
				reactionDict[reaction] = float(row["<Kd> (uM)"]) / 10**6


class QuoteDialect(csv.excel_tab):
	quotechar = "'"
	quoting = csv.QUOTE_NONNUMERIC

with io.open(REACTION_FILE, "rb") as infile:
	with io.open(OUTPUT_FILE, "wb") as outfile:
		reader = tsv.dict_reader(infile)
		fieldnames = list(reader.fieldnames or [])
		fieldnames.append("original reverse rate")
		writer = tsv.dict_writer(outfile, fieldnames=fieldnames, dialect=QuoteDialect)
		writer.writeheader()

		for row1 in reader:  # type: Dict[str, Any]
			row1["dir"] = float(row1["dir"])
			row1["forward rate"] = float(row1["forward rate"])
			row1["reverse rate"] = float(row1["reverse rate"])
			row1["stoichiometry"] = np.array(row1["stoichiometry"])
			if row1["id"] in reactionDict:
				row1["reverse rate"] = reactionDict[row1["id"]]
			row1["original reverse rate"] = row1["reverse rate"]

			writer.writerow(row1)
