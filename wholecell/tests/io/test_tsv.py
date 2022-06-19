"""Unit test for the tsv module."""
from __future__ import absolute_import, division, print_function

from io import BytesIO
from six.moves import zip_longest
import unittest

from wholecell.io import tsv


FIELD_NAMES = ['id', 'ourLocation', u'\u20ac:xyz', 'mass (units.g)']
INPUT_ROWS = [
	b'id\tourLocation\t\xE2\x82\xAC:xyz\tmass (units.g)',
	b'G6660-MONOMER\t[x]\tLocation information from Lopez Campistrous 2005.\t98.6',
	b'2.71828\t[c]\tLocation from \xe2\x8a\x972011.\t12']
FIELD_NAMES_WITH_PRIVATE_FIELD = ['id', 'ourLocation', u'_\u20ac:xyz', 'mass (units.g)']
INPUT_ROWS_WITH_PRIVATE_FIELD = [
	b'id\tourLocation\t_\xE2\x82\xAC:xyz\tmass (units.g)',
	b'G6660-MONOMER\t[x]\tLocation information from Lopez Campistrous 2005.\t98.6',
	b'2.71828\t[c]\tLocation from \xe2\x8a\x972011.\t12']

def _remove_private_fields(fieldnames):
	return [fieldname for fieldname in fieldnames if not fieldname.startswith('_')]

def _expected_row(index, private_field=False):
	values = INPUT_ROWS[index].decode('utf-8').split('\t')
	row_dict = {k: v for k, v in zip_longest(FIELD_NAMES, values, fillvalue=404)}

	if private_field:
		row_dict = {
			field_name: row_dict[field_name]
			for field_name in _remove_private_fields(FIELD_NAMES_WITH_PRIVATE_FIELD)
			}

	return row_dict


class Test_Tsv(unittest.TestCase):
	def test_reader(self):
		byte_stream = BytesIO(b'\n'.join(INPUT_ROWS))
		reader = tsv.reader(byte_stream)

		for row, expected in zip_longest(reader, INPUT_ROWS, fillvalue=404):
			assert row == expected.decode('utf-8').split('\t')

	def test_writer(self):
		byte_stream = BytesIO()
		writer = tsv.writer(byte_stream)
		for row in INPUT_ROWS:
			writer.writerow(row.decode('utf-8').split('\t'))

		data = byte_stream.getvalue()
		assert data == b'\r\n'.join(INPUT_ROWS + [b''])

	def test_dict_reader(self):
		byte_stream = BytesIO(b'\n'.join(INPUT_ROWS))
		reader = tsv.dict_reader(byte_stream)

		row_dict = next(reader)
		assert row_dict == _expected_row(1)

		row_dict = next(reader)
		assert row_dict == _expected_row(2)

		with self.assertRaises(StopIteration):
			next(reader)

		assert reader.fieldnames == FIELD_NAMES

		reader.fieldnames = FIELD_NAMES[1:]
		assert reader.fieldnames == FIELD_NAMES[1:]

	def test_dict_reader_with_private_field(self):
		byte_stream = BytesIO(b'\n'.join(INPUT_ROWS_WITH_PRIVATE_FIELD))
		reader = tsv.dict_reader(byte_stream)

		row_dict = next(reader)
		assert row_dict == _expected_row(1, private_field=True)

		row_dict = next(reader)
		assert row_dict == _expected_row(2, private_field=True)

		with self.assertRaises(StopIteration):
			next(reader)

		assert set(reader.fieldnames) == set(
			_remove_private_fields(FIELD_NAMES_WITH_PRIVATE_FIELD))

		reader.fieldnames = FIELD_NAMES_WITH_PRIVATE_FIELD[1:]
		assert set(reader.fieldnames) == set(
			_remove_private_fields(FIELD_NAMES_WITH_PRIVATE_FIELD[1:]))

	def test_dict_reader_with_initial_fieldnames(self):
		byte_stream = BytesIO(b'\n'.join(INPUT_ROWS) + b'\n')
		reader = tsv.dict_reader(byte_stream, fieldnames=FIELD_NAMES)

		row_dict = next(reader)
		assert row_dict == _expected_row(0)

		remaining_rows = list(reader)
		assert len(remaining_rows) == 2

	def test_dict_writer(self):
		field_names = FIELD_NAMES[2:]
		value2 = b'\xe2\x8a\x97xxx.'.decode('utf-8')

		byte_stream = BytesIO()
		writer = tsv.dict_writer(byte_stream, field_names, lineterminator='\n')

		writer.writeheader()
		writer.writerow({field_names[0]: 94022, field_names[1]: value2})
		writer.writerow({field_names[0]: '["c"]', field_names[1]: value2.upper()})

		data = byte_stream.getvalue()
		lines = data.split(b'\n')
		assert len(lines) == 4
		assert lines[0] == b'\xE2\x82\xAC:xyz\tmass (units.g)'
		assert lines[1] == b'94022\t' + value2.encode('utf-8')
		assert lines[2] == b'"[""c""]"\t' + value2.encode('utf-8').upper()
		assert lines[3] == b''
