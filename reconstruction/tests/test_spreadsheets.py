"""Unit test for the spreadsheets module."""
from __future__ import absolute_import, division, print_function

from io import BytesIO, TextIOWrapper
import os
import six
import unittest

from reconstruction.spreadsheets import JsonReader, JsonWriter, read_tsv
from wholecell.utils import units


FIELD_NAMES = ['id', 'ourLocation', u'\u20ac:xyz', 'mass (units.g)']
UNITLESS_FIELD_NAMES = FIELD_NAMES[:-1] + ['mass']
INPUT_DATA = b'''"id"\t"ourLocation"\t"\xE2\x82\xAC:xyz"\t"mass (units.g)"
"G6660-MONOMER"\t["c"]\t"Location information from Lopez Campistrous 2005."\t98.6
2.71828\t["c"]\t"Location from \xe2\x8a\x972011."\t12
'''


class Test_Spreadsheets(unittest.TestCase):
	def test_json_reader(self):
		byte_stream = BytesIO(INPUT_DATA)
		read_stream = byte_stream if six.PY2 else TextIOWrapper(byte_stream)
		reader = JsonReader(read_stream)
		l = list(reader)
		assert len(l) == 2
		assert l[0] == {
			'id': 'G6660-MONOMER',
			'ourLocation': ["c"],
			u'\u20ac:xyz': 'Location information from Lopez Campistrous 2005.',
			'mass': 98.6 * units.g}
		assert l[1] == {
			'id': 2.71828,
			'ourLocation': ["c"],
			u'\u20ac:xyz': b'Location from \xe2\x8a\x972011.'.decode('utf-8'),
			'mass': 12 * units.g}
		assert set(l[0].keys()) == set(UNITLESS_FIELD_NAMES)

	def test_read_tsv(self):
		filename = os.path.join('validation', 'ecoli', 'flat', 'schmidt2015_javier_table.tsv')
		entries = read_tsv(filename)
		assert b'Chemostat \xc2\xb5=0.20'.decode('utf-8') in entries[0]

	def test_json_writer(self):
		byte_stream = BytesIO()
		write_stream = byte_stream if six.PY2 else TextIOWrapper(byte_stream)
		field_names = UNITLESS_FIELD_NAMES  # JsonWriter doesn't support units!
		writer = JsonWriter(write_stream, field_names)
		writer.writeheader()

		row = {f: [f.upper()] for f in field_names}
		row['id'] = 3.14159
		row['mass'] = 1.23
		writer.writerow(row)

		data = byte_stream.getvalue()
		assert data == (
			b'"id"\t"ourLocation"\t"\xe2\x82\xac:xyz"\t"mass"\n'
			b'3.14159\t["OURLOCATION"]\t["\xe2\x82\xac:XYZ"]\t1.23\n')

	# TODO(jerry): Test dict and ndarray values.
	# TODO(jerry): Test tsv_writer() by writing and reading back a file.
