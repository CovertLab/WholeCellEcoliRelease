
import os
import re
import csv
import sys

import numpy as np

from reconstruction.ecoli.knowledge_base import KnowledgeBaseEcoli
from reconstruction.spreadsheets import JsonWriter

OUTPUT_DIR = sys.argv[1]

DIALECT = csv.excel_tab

kb = KnowledgeBaseEcoli(False)

drymass = kb._cellDryMassCompositionData

fraction_names = {
	"Glycogen":"glycogen",
	"InorganicIon":"ion",
	"LPS":"LPS",
	"Lipid":"lipid",
	"Murein":"murein",
	"SolublePool":"soluble"
	}

keys = drymass.dtype.names
with open(os.path.join(OUTPUT_DIR, "dryMassComposition.tsv"), "w") as outfile:
	writer = JsonWriter(outfile, keys, dialect = DIALECT)

	writer.writeheader()

	for entry in drymass:
		writer.writerow({
			key:entry[key] for key in keys
			})

for kb_name, file_name in fraction_names.viewitems():
	array = getattr(kb, "_cell{}FractionData".format(kb_name))
	keys = array.dtype.names

	with open(os.path.join(OUTPUT_DIR, "massFractions", "{}Fractions.tsv".format(file_name)), "w") as outfile:
		writer = JsonWriter(outfile, keys, dialect = DIALECT)

		writer.writeheader()

		for entry in array:
			writer.writerow({
				key:entry[key] for key in keys
				})
