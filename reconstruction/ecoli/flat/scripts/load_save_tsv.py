from __future__ import absolute_import, division, print_function

from reconstruction.spreadsheets import read_tsv, JsonWriter


def load_tsv(file_name):
	return read_tsv(file_name)

def write_tsv(list_of_dicts, file_name):
	with open(file_name, "w") as outfile:
		fieldnames = list_of_dicts[0].keys()
		writer = JsonWriter(outfile, fieldnames)
		writer.writeheader()
		writer.writerows(list_of_dicts)
