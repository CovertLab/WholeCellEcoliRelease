"""
CSV reader and writer that handle Python 2/3 compatibility and default to TAB
delimiters.
"""

from __future__ import absolute_import, division, print_function

import csv
from io import TextIOWrapper
from typing import Any, cast, IO, Iterable, List, Optional, Sequence, Text, Type, Union

import six


DIALECT = Union[str, csv.Dialect, Type[csv.Dialect]]


class reader(object):
	def __init__(self, csvfile, dialect='excel', delimiter='\t', **fmtparams):
		# type: (IO[bytes], DIALECT, str, **Any) -> None
		"""Open a csv reader(), handling Python 2/3 Unicode I/O compatibility,
		and defaulting to TAB delimiters.

		REQUIRES: `csvfile` must be a buffered byte reader, e.g. from
		io.open(filename, 'rb') or io.BytesIO(buffer).

		This does Unicode in Python 2 by decoding the values read by csv.reader,
		and in Python 3 by constructing the csv.reader with a TextIO.
		"""
		if six.PY2:
			self.input_file = csvfile
		else:
			self.input_file = TextIOWrapper(csvfile, encoding='utf-8', newline='')
		self.reader = csv.reader(
			self.input_file, dialect=dialect, delimiter=delimiter, **fmtparams)

	def __iter__(self):
		return self

	def __next__(self):
		# type: () -> List[Text]
		row = next(self.reader)
		if six.PY2:
			row_ = [value.decode('utf-8') for value in row]
			return row_
		return row

	next = __next__

	@property
	def dialect(self):
		# type: () -> DIALECT
		return self.reader.dialect

	@property
	def line_num(self):
		# type: () -> int
		return self.reader.line_num


class writer(object):
	def __init__(self, csvfile, dialect='excel', delimiter='\t', **fmtparams):
		# type: (IO[bytes], DIALECT, str, **Any) -> None
		"""Open a csv writer(), handling Python 2/3 Unicode I/O compatibility,
		and defaulting to TAB delimiters.

		REQUIRES: `csvfile` must be a buffered byte writer, e.g. from
		io.open(filename, 'wb') or io.BytesIO(buffer).

		This does Unicode in Python 2 by encoding the values passed to csv.writer,
		and in Python 3 by constructing the csv.writer with a TextIO.
		"""
		if six.PY2:
			self.output_file = csvfile
		else:
			self.output_file = TextIOWrapper(
				csvfile, encoding='utf-8', newline='', line_buffering=True)
		self.writer = csv.writer(
			self.output_file, dialect=dialect, delimiter=delimiter, **fmtparams)

	def writerow(self, row):
		# type: (Sequence[Any]) -> None
		def stringify(s):
			return s if isinstance(s, six.string_types) else Text(s)

		row_ = [stringify(value).encode('utf-8') for value in row] if six.PY2 else row
		self.writer.writerow(row_)

	def writerows(self, rows):
		# type: (Iterable[Sequence[Any]]) -> None
		for row in rows:
			self.writerow(row)

	@property
	def dialect(self):
		# type: () -> DIALECT
		return self.writer.dialect


def dict_reader(f, fieldnames=None, **kwargs):
	# type: (IO[bytes], Optional[List[str]], **Any) -> csv.DictReader
	"""Open a csv DictReader(), handling Python 2/3 Unicode I/O compatibility (by
	replacing its csv reader with a tsv.reader) and defaulting to TAB delimiters.

	REQUIRES: `csvfile` must be a buffered byte reader, e.g. from
	io.open(filename, 'rb') or io.BytesIO(buffer).
	"""
	tsv_reader = reader(f, **kwargs)
	tsv_dict_reader = csv.DictReader(
		tsv_reader.input_file, fieldnames=fieldnames, **kwargs)
	tsv_dict_reader.reader = cast(Any, tsv_reader)

	return tsv_dict_reader

def dict_writer(f, fieldnames, dialect='excel', **kwargs):
	# type: (IO[bytes], Iterable[str], DIALECT, **Any) -> csv.DictWriter
	"""Open a csv DictWriter(), handling Python 2/3 Unicode I/O compatibility (by
	replacing its csv writer with a tsv.writer) and defaulting to TAB delimiters.

	REQUIRES: `csvfile` must be a buffered byte writer, e.g. from
	io.open(filename, 'wb') or io.BytesIO(buffer).
	"""
	tsv_writer = writer(f, dialect=dialect, **kwargs)
	tsv_dict_writer = csv.DictWriter(
		tsv_writer.output_file, fieldnames, dialect=dialect, **kwargs)
	tsv_dict_writer.writer = cast(Any, tsv_writer)
	# TODO(jerry): Call `tsv_dict_writer.writeheader()` for convenience?

	return tsv_dict_writer
