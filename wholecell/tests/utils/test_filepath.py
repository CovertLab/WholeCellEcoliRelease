"""
Test the filepath utility.
"""

from __future__ import absolute_import, division, print_function

import os
import shutil
import tempfile
import unittest

import pytest

from wholecell.utils import filepath


class Test_filepath(unittest.TestCase):

	def setUp(self):
		# Create a temporary directory
		self.test_dir = tempfile.mkdtemp()

	def tearDown(self):
		shutil.rmtree(self.test_dir)

	def test_makedirs(self):
		"""Test makedirs(), write_file(), verify_file_exists(), and
		verify_dir_exists().
		"""
		directories = 'this/is/a/test'
		expected_path = os.path.join(self.test_dir, directories)
		self.assertFalse(os.path.exists(expected_path))

		# Test creating a directory path.
		with pytest.raises(IOError):
			filepath.verify_file_exists(expected_path)
		with pytest.raises(IOError):
			filepath.verify_dir_exists(expected_path)
		full_path = filepath.makedirs(self.test_dir, directories)
		self.assertEqual(full_path, expected_path)
		self.assertTrue(os.path.exists(expected_path))

		# Test that it's happy with an existing path.
		full_path2 = filepath.makedirs(self.test_dir, 'this', 'is', 'a/test')
		self.assertEqual(full_path2, expected_path)
		self.assertTrue(os.path.exists(expected_path))
		with pytest.raises(IOError):
			filepath.verify_file_exists(expected_path)
		filepath.verify_dir_exists(expected_path)

		# Test failure to create a directory path because a data file is there.
		filename = 'data'
		output_path = os.path.join(full_path, filename)
		with pytest.raises(IOError):
			filepath.verify_file_exists(expected_path)
		filepath.write_file(output_path, 'hi')
		filepath.verify_file_exists(output_path)
		with pytest.raises(OSError):
			filepath.makedirs(self.test_dir, directories, filename)

	def test_timestamp(self):
		"""Test timestamp() and TIMESTAMP_PATTERN."""
		timestamp = filepath.timestamp()
		self.assertTrue(timestamp.startswith('20'),
			'Timestamp "{}" starts with "20"'.format(timestamp))
		self.assertEqual(len(timestamp), 15, 'len(timestamp)')
		self.assertRegexpMatches(timestamp, filepath.TIMESTAMP_PATTERN)

	def test_json_files(self):
		"""Test read_json_file(), write_json_file."""
		expected = {'a': range(10), 'b': 14.5, 'c': None}
		output_path = os.path.join(self.test_dir, 'abc.json')
		filepath.write_json_file(output_path, expected)
		actual = filepath.read_json_file(output_path)
		self.assertEqual(expected, actual)

	def test_iter_variants(self):
		results = list(filepath.iter_variants('XYZ', 123, 125))
		self.assertEqual(
			[(123, 'XYZ_000123'), (124, 'XYZ_000124'), (125, 'XYZ_000125')],
			results)


if __name__ == '__main__':
	unittest.main()
