"""Tools for working with dependency graphs"""


from typing import List, Dict, Iterable


class InvalidDependencyGraphError(Exception):
	"""Exception for invalid dependency graphs"""
	pass


class DependencyGraph(object):
	"""Represents a dependency graph"""

	#: Mapping from node to that node's dependents
	dependencies = {}  # type: Dict[str, List[str]]

	def __init__(self):
		"""Initialize dependencies to empty dictionary"""
		self.dependencies = dict()

	def add_nodes(self, nodes):
		# type: (List[str]) -> None
		"""Add nodes with no dependencies

		Arguments:
			nodes: Nodes to add
		"""
		# type: (Iterable[str]) -> None
		for node in nodes:
			self.dependencies[node] = []

	def add_dep_relation(self, a, b):
		# type: (str, str) -> None
		"""Add an edge such that a depends on b"""
		if a in self.dependencies:
			self.dependencies[a].append(b)
		else:
			self.dependencies[a] = [b]
		if b not in self.dependencies:
			self.dependencies[b] = []

	def get_topological_ordering(self):
		# type: () -> List[str]
		"""Get a topological ordering of the nodes

		Returns:
			List of denendency names such that the dependencies can be
			loaded in the order in which they appear in the list without
			violating dependency relationships.

		Raises:
			InvalidDependencyGraphError: If the graph contains a cycle
		"""
		explored = {name: False for name in self.dependencies}
		rev_ordering = []  # type: List[str]
		for node in self.dependencies:
			if not explored[node]:
				self._topo_sort_dfs(node, explored, rev_ordering)
		return rev_ordering

	def _topo_sort_dfs(self, node, explored, rev_ordering):
		# type: (str, Dict[str, bool], List[str]) -> None
		explored[node] = True
		for dependency in self.dependencies[node]:
			if not explored[dependency]:
				self._topo_sort_dfs(dependency, explored, rev_ordering)
			else:
				if dependency not in rev_ordering:
					# We have reached a node that we have already
					# visited, but that we have not added to the
					# ordering. Thus a call to dependency is still on
					# the stack. This means there exists a path from
					# dependency to node, so there is a cycle.
					raise InvalidDependencyGraphError(
						"Dependency graphs must not contain cycles")
		rev_ordering.append(node)
