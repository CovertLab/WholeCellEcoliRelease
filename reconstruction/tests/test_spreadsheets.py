from __future__ import absolute_import, division, print_function

import csv
from io import BytesIO, TextIOWrapper
import unittest

from reconstruction.spreadsheets import JsonReader, JsonWriter


FIELD_NAMES = ['id', 'ourLocation', 'comments', 'ecocycLocations']
INPUT_DATA = b'''"id"\t"ourLocation"\t"comments"\t"ecocycLocations"
"G6660-MONOMER"\t["c"]\t"Location information from Lopez Campistrous 2005."\t["z"]
"FHUC-MONOMER"\t["c"]\t"Location information from Han 2011."\t["i"]
'''


CSV_DIALECT = csv.excel_tab


class Test_Spreadsheets(unittest.TestCase):
	def test_reader(self):
		byte_stream = BytesIO(INPUT_DATA)
		text_stream = TextIOWrapper(byte_stream)
		reader = JsonReader(text_stream, dialect=CSV_DIALECT)
		l = list(reader)
		assert len(l) == 2
		assert l[0] == {
			'id': 'G6660-MONOMER',
			'ourLocation': ["c"],
			'comments': 'Location information from Lopez Campistrous 2005.',
			'ecocycLocations': ['z']}
		assert l[1] == {
			'id': 'FHUC-MONOMER',
			'ourLocation': ["c"],
			'comments': 'Location information from Han 2011.',
			'ecocycLocations': ['i']}
		assert set(l[0].keys()) == set(FIELD_NAMES)

	def test_writer(self):
		byte_stream = BytesIO(b'')
		# TODO(jerry): JsonWriter might need to write unicode to a text stream
		#  for Python 3 compatibility, but it writes str in Python 2.
		#  Is that a bug in JsonWriter or fighting with the Python 2 & 3
		#  libraries? E.g. Python 3 doesn't have `StringIO.StringIO`.
		# text_stream = TextIOWrapper(byte_stream)
		writer = JsonWriter(byte_stream, FIELD_NAMES, dialect=CSV_DIALECT)
		writer.writeheader()

		writer.writerow({f: [f.upper()] for f in FIELD_NAMES})

		data = byte_stream.getvalue()
		assert data == (
			b'"id"\t"ourLocation"\t"comments"\t"ecocycLocations"\n'
			'["ID"]\t["OURLOCATION"]\t["COMMENTS"]\t["ECOCYCLOCATIONS"]\n')
