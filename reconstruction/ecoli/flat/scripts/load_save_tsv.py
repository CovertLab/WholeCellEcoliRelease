from __future__ import absolute_import, division, print_function

from reconstruction.spreadsheets import read_tsv, tsv_writer


def load_tsv(file_name):
	return read_tsv(file_name)

def write_tsv(list_of_dicts, file_name):
	fieldnames = list(list_of_dicts[0].keys())
	with tsv_writer(file_name, fieldnames) as writer:
		writer.writerows(list_of_dicts)
