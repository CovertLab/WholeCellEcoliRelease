"""Unit test for the spreadsheets module."""
from __future__ import absolute_import, division, print_function

from io import BytesIO, TextIOWrapper
import os
import six
import unittest

from reconstruction.spreadsheets import JsonReader, JsonWriter, read_tsv, tsv_reader
from wholecell.utils import units


JAVIER_TABLE = os.path.join(
	'validation', 'ecoli', 'flat', 'schmidt2015_javier_table.tsv')
CHEMOSTAT_20 = b'Chemostat \xc2\xb5=0.20'.decode('utf-8')

FIELD_NAMES = ['id', 'ourLocation', u'\u20ac:xyz', 'mass (units.g)']
UNITLESS_FIELD_NAMES = FIELD_NAMES[:-1] + ['mass']
INPUT_DATA = b'''"id"\t"ourLocation"\t"\xE2\x82\xAC:xyz"\t"mass (units.g)"
"G6660-MONOMER"\t["c"]\t"Location information from Lopez Campistrous 2005."\t98.6
2.71828\t["c"]\t"Location from \xe2\x8a\x972011."\t12
'''
INPUT_DATA_WITH_PRIVATE_FIELD = b'''"id"\t"ourLocation"\t"_\xE2\x82\xAC:xyz"\t"mass (units.g)"
"G6660-MONOMER"\t["c"]\t"Location information from Lopez Campistrous 2005."\t98.6
2.71828\t["c"]\t"Location from \xe2\x8a\x972011."\t12
'''
NONPRIVATE_FIELD_NAMES = ['id', 'ourLocation', 'mass']


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

	def test_json_reader_with_private_field(self):
		byte_stream = BytesIO(INPUT_DATA_WITH_PRIVATE_FIELD)
		read_stream = byte_stream if six.PY2 else TextIOWrapper(byte_stream)
		reader = JsonReader(read_stream)
		l = list(reader)
		assert len(l) == 2
		assert l[0] == {
			'id': 'G6660-MONOMER',
			'ourLocation': ["c"],
			'mass': 98.6 * units.g}
		assert l[1] == {
			'id': 2.71828,
			'ourLocation': ["c"],
			'mass': 12 * units.g}
		assert set(l[0].keys()) == set(NONPRIVATE_FIELD_NAMES)

	def test_tsv_reader(self):
		with tsv_reader(JAVIER_TABLE) as reader:
			fieldnames = reader.fieldnames
			assert fieldnames[0] == 'EcoCycID'
			assert fieldnames[-1] == 'Fructose'
			row1 = next(reader)

		assert row1['EcoCycID'] == 'EG10001'
		assert row1[CHEMOSTAT_20] == 9
		assert CHEMOSTAT_20 in row1

		if six.PY2:
			# reader.fieldnames contains UTF-8 `bytes`. rows keys are unicode strings.
			assert CHEMOSTAT_20.encode('utf-8') in fieldnames
		else:
			assert CHEMOSTAT_20 in fieldnames

	def test_read_tsv(self):
		entries = read_tsv(JAVIER_TABLE)
		assert CHEMOSTAT_20 in entries[0]

	def test_json_writer(self):
		byte_stream = BytesIO()
		write_stream = byte_stream if six.PY2 else TextIOWrapper(
			byte_stream, line_buffering=True)
		key2 = u'key \u20ac.'
		field_names = ['key1', key2, '33.3']
		writer = JsonWriter(write_stream, field_names)
		writer.writeheader()

		writer.writerow({'key1': 'value1', key2: '', '33.3': 33.4})
		writer.writerow({'key1': [1.1, '11'], key2: None, '33.3': {u'\u2297!': [33]}})

		data = byte_stream.getvalue()
		lines = data.split(b'\n')
		assert len(lines) == 4
		assert lines[0] == b'"key1"\t"key \xE2\x82\xAC."\t"33.3"'
		assert lines[1] == b'"value1"\t""\t33.4'
		assert lines[2] == b'[1.1, "11"]\tnull\t{"\xe2\x8a\x97!": [33]}'
		assert lines[3] == b''

	# TODO(jerry): Test dict and ndarray values.
	# TODO(jerry): Test tsv_writer() by writing and reading back a file.
