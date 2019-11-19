#!/usr/bin/env python2

"""Tests for behavior metrics"""

from __future__ import absolute_import, division, print_function
from typing import Dict
import unittest

import mock
import numpy as np

from wholecell.io.tablereader import TableReader
from prototypes.behavior_metrics.behavior_metrics import BehaviorMetrics
from prototypes.behavior_metrics.dependency_graph import (
	InvalidDependencyGraphError
)


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

	def test_order_operations(self):
		# Graph where an edge A -> B means that B depends on A
		# A -> B
		#       \
		#        +-> D
		#       /
		#      C
		config = {
			"A": {},
			"B": {"args": ["A"]},
			"C": {},
			"D": {"args": ["B", "C"]}
		}
		order = BehaviorMetrics.order_operations(config)
		indexed_order = {source: i for i, source in enumerate(order)}
		self._assertComesBefore(indexed_order, "A", "B")
		self._assertComesBefore(indexed_order, "B", "D")
		self._assertComesBefore(indexed_order, "C", "D")
		self.assertEqual(4, len(order))

	def test_order_ops_cycle(self):
		# Graph where an edge A -> B means that B depends on A
		# A --> B
		# ^     |
		# |     |
		# +- C <+
		config = {
			"A": {"args": ["C"]},
			"B": {"args": ["A"]},
			"C": {"args": ["B"]},
		}
		with self.assertRaisesRegexp(InvalidDependencyGraphError, "cycle"):
			BehaviorMetrics.order_operations(config)

	def _assertComesBefore(self, indexed_order, earlier, later):
		# type: (Dict[str, int], str, str) -> None
		self.assertLess(
			indexed_order[earlier], indexed_order[later],
			"{} does not come before {}".format(earlier, later)
		)


if __name__ == "__main__":
	unittest.main()
