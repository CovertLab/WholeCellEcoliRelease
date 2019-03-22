#!/usr/bin/env python

# Generates a tsv output file to replace equilibriumReactions.tsv with updated Kd values
# Requires 2 files: validation/ecoli/flat/conformationTF.tsv and reconstruction/ecoli/flat/equilibriumReactions.tsv
# Outputs TFoutput.tsv to directory that script is run from
# Need to replace ' with " in file after script runs

from __future__ import division

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import wholecell.utils.constants
from wholecell.io.tablereader import TableReader

import csv
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

CSV_DIALECT = csv.excel_tab
DATA_FILE = os.path.join("validation","ecoli","flat","conformationTF.tsv")
REACTION_FILE = os.path.join("reconstruction","ecoli","flat","equilibriumReactions.tsv")
OUTPUT_FILE = "TFoutput.tsv"

reactionDict = {}
with open(DATA_FILE, "rU") as csvfile:
	reader = csv.DictReader(csvfile, dialect = CSV_DIALECT)
	for row in reader:
		if row["<Kd> (uM)"] != '' and row["<Kd> (uM)"] != '?':
			for reaction in row["EcoCyc ID reaction (metabolite vs. TF)"].split(", "):
				reactionDict[reaction] = float(row["<Kd> (uM)"]) / 10**6

with open(REACTION_FILE, "rU") as infile:
	with open(OUTPUT_FILE, "w") as outfile:
		reader = csv.DictReader(infile, dialect = CSV_DIALECT)
		quoteDialect = CSV_DIALECT
		quoteDialect.quotechar = "'"
		quoteDialect.quoting = csv.QUOTE_NONNUMERIC
		writer = csv.DictWriter(outfile, fieldnames = reader.fieldnames, dialect = quoteDialect)
		writer.fieldnames.append("original reverse rate")
		writer.writeheader()
		for row in reader:
			row["dir"] = float(row["dir"])
			row["forward rate"] = float(row["forward rate"])
			row["reverse rate"] = float(row["reverse rate"])
			row["stoichiometry"] = np.array(row["stoichiometry"])
			if row["id"] in reactionDict.keys():
				row["reverse rate"] = reactionDict[row["id"]]
			row["original reverse rate"] = row["reverse rate"]

			writer.writerow(row)
