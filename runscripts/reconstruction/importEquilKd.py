#!/usr/bin/env python

# Generates a tsv output file to replace equilibriumReactions.tsv with updated Kd values
# Requires 2 files: validation/ecoli/flat/conformationTF.tsv and reconstruction/ecoli/flat/equilibriumReactions.tsv
# Outputs TFoutput.tsv to directory that script is run from
# Need to replace ' with " in file after script runs

from __future__ import absolute_import, division, print_function

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")

import csv
from typing import Any, Dict


CSV_DIALECT = csv.excel_tab
DATA_FILE = os.path.join("validation","ecoli","flat","conformationTF.tsv")
REACTION_FILE = os.path.join("reconstruction","ecoli","flat","equilibriumReactions.tsv")
OUTPUT_FILE = "TFoutput.tsv"

reactionDict = {}
with open(DATA_FILE, "rU") as csvfile:
	reader = csv.DictReader(csvfile, dialect = CSV_DIALECT)
	for row in reader:  # type: Dict[str, Any]
		if row["<Kd> (uM)"] != '' and row["<Kd> (uM)"] != '?':
			for reaction in row["EcoCyc ID reaction (metabolite vs. TF)"].split(", "):
				reactionDict[reaction] = float(row["<Kd> (uM)"]) / 10**6

with open(REACTION_FILE, "rU") as infile:
	with open(OUTPUT_FILE, "w") as outfile:
		reader = csv.DictReader(infile, dialect = CSV_DIALECT)
		quoteDialect = CSV_DIALECT
		quoteDialect.quotechar = "'"
		quoteDialect.quoting = csv.QUOTE_NONNUMERIC

		fieldnames = list(reader.fieldnames or [])
		fieldnames.append("original reverse rate")
		writer = csv.DictWriter(outfile, fieldnames=fieldnames, dialect=quoteDialect)
		writer.writeheader()

		for row1 in reader:  # type: Dict[str, Any]
			row1["dir"] = float(row1["dir"])
			row1["forward rate"] = float(row1["forward rate"])
			row1["reverse rate"] = float(row1["reverse rate"])
			row1["stoichiometry"] = np.array(row1["stoichiometry"])
			if row1["id"] in reactionDict.keys():
				row1["reverse rate"] = reactionDict[row1["id"]]
			row1["original reverse rate"] = row1["reverse rate"]

			writer.writerow(row1)
