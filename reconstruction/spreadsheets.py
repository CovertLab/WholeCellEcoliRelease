"""
Subclasses of DictWriter and DictReader that parse plaintext as JSON strings,
allowing for basic type parsing and fields that are dictionaries or lists. The
reader also supports units and comment lines.
"""

from __future__ import absolute_import, division, print_function

from contextlib import contextmanager
import csv
import io
from itertools import filterfalse
import json
import re
import numpy as np
from typing import Any, Dict, Iterator, List, Sequence

from wholecell.utils import units


CSV_DIALECT = csv.excel_tab


def comment_line(line):
	# type: (str) -> bool
	return line.lstrip().startswith('#')


# TODO(jerry): Rename tsv_reader() et al to tsv_json_reader()?
# TODO(jerry): Implementing this on wholecell.io.tsv.dict_reader/dict_writer
#  would simplify it, but not as much after we delete the PY2 code.

@contextmanager
def tsv_reader(filename):
	# type: (str) -> Iterator[JsonReader]
	"""A context manager that opens a TSV JsonReader on the given filename with
	an input filter to skip comment lines.
	"""
	# ########################################################################
	# NOTE: Python 3 csv requires opening the file as text 'r' with the right
	# character encoding while Python 2 csv requires opening it as bytes 'rb'
	# then decoding from UTF-8 after csv reads it but before DictReader or at
	# least before json.loads(). The file can be in UTF-8 or its ASCII subset.
	# ########################################################################

	with io.open(filename, mode='r', encoding='utf-8', newline='') as fh:
		reader = JsonReader(filterfalse(comment_line, fh), dialect=CSV_DIALECT)
		yield reader


def read_tsv(filename):
	# type: (str) -> Sequence[Dict[str, Any]]
	"""Read an entire .tsv file using JsonReader and skip comment lines."""
	with tsv_reader(filename) as reader:
		return list(reader)


@contextmanager
def tsv_writer(filename, fieldnames):
	# type: (str, Sequence[str]) -> Iterator[JsonWriter]
	"""A context manager that opens a TSV JsonWriter on the given filename.
	Just call its writerow() and writerows() methods.
	"""
	# ########################################################################
	# NOTE: Python 3 csv requires opening the file as text 'w' with the right
	# character encoding while Python 2 csv requires opening it as bytes 'wb'
	# and encoding the data to UTF-8 before csv writes it.
	# ########################################################################

	fieldnames = list(fieldnames)
	with io.open(filename, mode='w', encoding='utf-8', newline='') as fh:
		writer = JsonWriter(fh, fieldnames, dialect=CSV_DIALECT)
		writer.writeheader()

		yield writer

		fh.flush()


def array_to_list(value):
	if isinstance(value, np.ndarray):
		value = value.tolist()

	return value


class JsonWriter(csv.DictWriter):
	def __init__(self, *args, **kwargs):
		"""Writer for a .tsv file to be read by JsonReader. This writes a
		header with quotes and dict rows in TSV format, JSON encoding, and
		UTF-8 encoding.

		NOTE: The caller needs to remove units from the dict values and add
		them to the fieldnames. JsonWriter does not handle units.

		The first argument should be a file-like writer open in text mode.

		By default, dialect=CSV_DIALECT, which is excel-tab.
		"""
		kwargs.setdefault('dialect', CSV_DIALECT)
		super(JsonWriter, self).__init__(
			quotechar = "'", quoting = csv.QUOTE_MINIMAL, lineterminator="\n", *args, **kwargs
			)

	def writeheader(self):
		# Bypass DictWriter's writeheader() and _dict_to_list(). [Consider
		# reimplementing on csv.writer rather than subclassing DictWriter.]
		header = [u'"{}"'.format(name) for name in self.fieldnames]
		self.writer.writerow(header)

	def _dict_to_list(self, rowdict):
		rowdict_ = {
			key: json.dumps(array_to_list(value), ensure_ascii=False)
			for key, value in rowdict.items()}
		# noinspection PyUnresolvedReferences
		return super(JsonWriter, self)._dict_to_list(rowdict_)


class JsonReader(object):
	def __init__(self, *args, **kwargs):
		"""
		Reader for a .tsv file that supports units and json-coded values.
		Units are denoted with a fieldname in the format 'name (units)' e.g.
		"flux standard deviation (units.mmol / units.g / units.h)". Fields
		whose names start with an underscore are removed from self._fieldnames,
		and discarded from each row during iteration.

		The first argument should be a file-like reader open in text mode.

		By default, dialect=CSV_DIALECT, which is excel-tab.
		"""
		kwargs.setdefault('dialect', CSV_DIALECT)
		self.tsv_dict_reader = csv.DictReader(
			quotechar = "'", quoting = csv.QUOTE_MINIMAL, *args, **kwargs
			)

		fieldnames = self.tsv_dict_reader.fieldnames

		# Strip extra quotes from the field names
		fieldnames = [fieldname.strip('"') for fieldname in fieldnames]
		self.tsv_dict_reader.fieldnames = fieldnames

		# Discard private field names that begin with underscore and empty field names
		self._fieldnames = [
			fieldname for fieldname in fieldnames
			if not fieldname.startswith('_') and fieldname != '']

	def __iter__(self):
		return self

	def __next__(self):
		# type: () -> Dict[str, Any]
		return self._decode_row(self.tsv_dict_reader.__next__())

	@property
	def fieldnames(self):
		# type: () -> List[str]
		return self._fieldnames

	def _decode_row(self, row_dict):
		# type: (Dict[str, str]) -> Dict[str, Any]
		"""Decode a DictReader row.

		NOTE: Each returned row contains unicode/str keys and values.
		"""
		attributeDict = {}  # type: Dict[str, Any]

		for fieldname in self._fieldnames:
			raw_value = row_dict[fieldname]
			key = fieldname

			try:
				value = json.loads(raw_value) if raw_value else ""

			except (ValueError, TypeError) as e:
				repr(e)  # TODO(jerry): Why call repr() and discard the result?
				raise ValueError("failed to parse json string:{}".format(raw_value))

			match = re.search(r'(.*?) \((.*?)\)', key)
			if match:
				# Entry has units so apply the units to the value and strip
				# them from the key.
				_ = units  # don't warn about `units`; it's imported for eval()
				attribute = match.group(1)
				value_units = eval(match.group(2))

				# Units do not work with empty values
				if value != '':
					if isinstance(value, dict):
						# Apply units to each dictionary value
						for k, v in value.items():
							value[k] = value_units * v
							value[k].normalize()
					else:
						# Units do not work with lists so need to convert to ndarray
						if isinstance(value, list):
							value_with_units = value_units * np.array(value)
						else:
							value_with_units = value_units * value
						# Normalize to catch any unit issues now instead of later
						value = value_with_units.normalize()

				attributeDict[attribute] = value
			else:
				attributeDict[key] = value

		return attributeDict

	@property
	def dialect(self):
		# type: () -> CSV_DIALECT
		return self.tsv_dict_reader.dialect

	@property
	def line_num(self):
		# type: () -> int
		return self.tsv_dict_reader.line_num
