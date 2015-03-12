
import os
import re
import csv
import sys

import numpy as np

from reconstruction.ecoli.knowledge_base import KnowledgeBaseEcoli
from reconstruction.spreadsheets import JsonWriter

OUTPUT_DIR = sys.argv[1]
REPORT_FILENAME = os.path.join(OUTPUT_DIR, "generated_file_report.txt")
DIALECT = csv.excel_tab

kb = KnowledgeBaseEcoli(False)

report_file = sys.stdout #open(REPORT_FILENAME, "w")

candidates = []

for attr_name in dir(kb):

	action = None

	attr = getattr(kb, attr_name)

	if re.match("^_[a-zA-Z]\w+", attr_name) is None:
		action = "ignored because of pattern: "

	elif callable(attr):
		action = "ignored because callable: "

	elif not isinstance(attr, list):
		action = "ignored because not a list: "

	elif not isinstance(attr[0], dict):
		action = "ignored because values are not dicts: "

	else:
		action = "added as candidate: "
		candidates.append(attr_name)

	report_file.write(
		action + attr_name + "\n"
		)

for attr_name in candidates:
	try:
		attr = getattr(kb, attr_name)

		fieldnames = attr[0].keys()

		outfile = open(os.path.join(OUTPUT_DIR, attr_name.lstrip("_")) + ".tsv", "w")

		writer = JsonWriter(outfile, fieldnames, dialect = DIALECT)

		writer.writeheader()

		writer.writerows(attr)

	except Exception as e:
		# report_file.write("failed to write {} with exception {}\n".format(
		# 	attr_name,
		# 	e.message
		# 	))
		raise

	else:
		report_file.write("wrote {} to file\n".format(attr_name))

report_file.close()
