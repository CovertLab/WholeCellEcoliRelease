"""
Test the filepath utility.
"""

from __future__ import absolute_import
from __future__ import division

import nose.plugins.attrib as noseAttrib
import nose.tools
import os
import shutil
import tempfile
import unittest

from wholecell.utils import filepath

class Test_filepath(unittest.TestCase):

	def setUp(self):
		# Create a temporary directory
		self.test_dir = tempfile.mkdtemp()

	def tearDown(self):
		shutil.rmtree(self.test_dir)

	@noseAttrib.attr('smalltest', 'filepath')
	def test_makedirs(self):
		directories = 'this/is/a/test'
		expected_path = os.path.join(self.test_dir, directories)
		self.assertFalse(os.path.exists(expected_path))

		# Test creating a directory path.
		full_path = filepath.makedirs(self.test_dir, directories)
		self.assertEqual(full_path, expected_path)
		self.assertTrue(os.path.exists(expected_path))

		# Test that it's happy with an existing path.
		full_path2 = filepath.makedirs(self.test_dir, 'this', 'is', 'a/test')
		self.assertEqual(full_path2, expected_path)
		self.assertTrue(os.path.exists(expected_path))

		# Test failure to create a directory path because a data file is there.
		filename = 'data'
		with open(os.path.join(full_path, filename), 'w') as f:
			f.write('hi')
		with nose.tools.assert_raises(OSError):
			filepath.makedirs(self.test_dir, directories, filename)


if __name__ == '__main__':
	unittest.main()
