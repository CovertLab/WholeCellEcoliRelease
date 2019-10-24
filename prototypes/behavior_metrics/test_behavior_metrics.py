#!/usr/bin/env python3

import tempfile
import unittest

import mock
import numpy as np

from wholecell.io.tablereader import TableReader
from prototypes.behavior_metrics.behavior_metrics import BehaviorMetrics


# pragma pylint: disable=missing-docstring


class TestParseDataConfig(unittest.TestCase):

	maxDiff = None

	def setUp(self):
		self.mock_col = np.ndarray([1, 2, 3])
		self.mock_attr = 5
		self.mock_reader = mock.MagicMock(
			readColumn=mock.MagicMock(
				return_value=self.mock_col,
			),
			readAttribute=mock.MagicMock(
				return_value=self.mock_attr,
			),
		)
		self.mock_tablereader = mock.patch(
			"prototypes.behavior_metrics.behavior_metrics.TableReader",
			return_value=self.mock_reader
		).start()

	def test_load_column(self):
		config_str = {"A": {"table": "myTable", "column": "myCol"}}
		behavior_metrics = BehaviorMetrics("myConfPath", "mySimOutDir")
		data = behavior_metrics.load_data_from_config(config_str)
		self.mock_tablereader.assert_called_with("mySimOutDir/myTable")
		self.mock_reader.readColumn.assert_called_with("myCol")
		self.assertDictEqual({"A": self.mock_col}, data)

	def test_load_attribute(self):
		config_str = {"A": {"table": "myTable", "attribute": "myAttr"}}
		behavior_metrics = BehaviorMetrics("myConfPath", "mySimOutDir")
		data = behavior_metrics.load_data_from_config(config_str)
		self.mock_tablereader.assert_called_with("mySimOutDir/myTable")
		self.mock_reader.readAttribute.assert_called_with("myAttr")
		self.assertDictEqual({"A": self.mock_attr}, data)


if __name__ == "__main__":
	unittest.main()
