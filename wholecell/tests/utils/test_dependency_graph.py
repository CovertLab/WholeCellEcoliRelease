#!/usr/bin/env python3


"""Tests for Dependency Graph"""

from __future__ import absolute_import, division, print_function

import unittest

from wholecell.utils.dependency_graph import (
	DependencyGraph,
	InvalidDependencyGraphError,
)


class TestDependencyGraph(unittest.TestCase):

	maxDiff = None

	def test_graph_creation(self):
		# Graph:
		#     +--> B --+
		#     |        |
		# A --+        +--> D
		#     |        |
		#     +--> C --+
		dg = DependencyGraph()
		dg.add_dep_relation("B", "A")
		dg.add_dep_relation("C", "A")
		dg.add_dep_relation("D", "B")
		dg.add_dep_relation("D", "C")

		expected = {
			"D": ["B", "C"],
			"B": ["A"],
			"C": ["A"],
			"A": [],
		}
		self.assertDictEqual(expected, dg.dependencies)

	def test_topological_ordering_single_path(self):
		# Graph: A -> B -> C
		dg = DependencyGraph()
		dg.add_dep_relation("C", "B")
		dg.add_dep_relation("B", "A")
		ordering = dg.get_topological_ordering()
		self.assertEqual(["A", "B", "C"], ordering)

	def test_topological_ordering_multi_path(self):
		# Graph:
		#     +--> B --+
		#     |        |
		# A --+        +--> D
		#     |        |
		#     +--> C --+
		dg = DependencyGraph()
		dg.add_dep_relation("B", "A")
		dg.add_dep_relation("C", "A")
		dg.add_dep_relation("D", "B")
		dg.add_dep_relation("D", "C")

		ordering = dg.get_topological_ordering()

		self.assertIn(ordering, [
			["A", "B", "C", "D"],
			["A", "C", "B", "D"],
		])

	def test_topological_ordering_empty(self):
		dg = DependencyGraph()
		ordering = dg.get_topological_ordering()
		self.assertListEqual([], ordering)

	def test_topological_ordering_single_node(self):
		dg = DependencyGraph()
		dg.add_nodes(["A"])
		ordering = dg.get_topological_ordering()
		self.assertListEqual(["A"], ordering)

	def test_topological_ordering_non_hierarchical(self):
		# Graph:
		#     +--------+
		#     |        |
		# A --+        +--> D
		#     |        |
		#     +--> C --+
		dg = DependencyGraph()
		dg.add_dep_relation("D", "A")
		dg.add_dep_relation("D", "C")
		dg.add_dep_relation("C", "A")

		ordering = dg.get_topological_ordering()
		self.assertListEqual(["A", "C", "D"], ordering)

	def test_topological_ordering_cycle(self):
		dg = DependencyGraph()
		dg.add_dep_relation("A", "B")
		dg.add_dep_relation("B", "C")
		dg.add_dep_relation("C", "A")

		with self.assertRaisesRegexp(
			InvalidDependencyGraphError, "cycle"
		):
			dg.get_topological_ordering()

	def test_topological_ordering_unconnected(self):
		dg = DependencyGraph()
		nodes = ["A", "B", "C"]
		dg.add_nodes(nodes)
		ordering = dg.get_topological_ordering()
		self.assertSetEqual(set(nodes), set(ordering))
		self.assertEqual(len(nodes), len(ordering))


if __name__ == "__main__":
	unittest.main()
