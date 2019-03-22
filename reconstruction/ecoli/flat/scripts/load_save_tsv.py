import os
import csv
from reconstruction.spreadsheets import JsonReader, JsonWriter


def load_tsv(file_name):
	with open(file_name, 'rU') as csvfile:
		reader = JsonReader(csvfile, dialect = csv.excel_tab)
		return [row for row in reader]

def write_tsv(list_of_dicts, file_name):
	with open(file_name, "w") as outfile:
		fieldnames = list_of_dicts[0].keys()
		writer = JsonWriter(outfile, fieldnames, dialect = csv.excel_tab)
		writer.writeheader()
		writer.writerows(list_of_dicts)