from __future__ import absolute_import, division, print_function

import os
import csv
import sys

import six

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from reconstruction.spreadsheets import tsv_writer

OUTPUT_DIR = sys.argv[1]

DIALECT = csv.excel_tab

kb = KnowledgeBaseEcoli()

drymass = getattr(kb, '_cellDryMassCompositionData')

fraction_names = {
	"Glycogen":"glycogen",
	"InorganicIon":"ion",
	"LPS":"LPS",
	"Lipid":"lipid",
	"Murein":"murein",
	"SolublePool":"soluble"
	}

keys = drymass.dtype.names
filename = os.path.join(OUTPUT_DIR, "dry_mass_composition.tsv")
with tsv_writer(filename, keys) as writer:
	for entry in drymass:
		writer.writerow({
			key:entry[key] for key in keys
			})

for kb_name, file_prefix in six.viewitems(fraction_names):
	array = getattr(kb, "_cell{}FractionData".format(kb_name))
	keys = array.dtype.names

	filename = os.path.join(
		OUTPUT_DIR, "mass_fractions", "{}Fractions.tsv".format(file_prefix))
	with tsv_writer(filename, keys) as writer:
		for entry in array:
			writer.writerow({
				key:entry[key] for key in keys
				})
