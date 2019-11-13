"""Tools for working with dependency graphs"""


from typing import List, Dict, Iterable


class InvalidDependencyGraphError(Exception):
	"""Exception for invalid dependency graphs"""
	pass


class ExplorationStatus(object):
	#: Node has not yet been visited.
	UNEXPLORED = 1
	#: Node has been visited but not added to list. This means that we
	#: have not yet finished with the node.
	EXPLORING = 2
	#: Node has been visited and added to list.
	EXPLORED = 3


class DependencyGraph(object):
	"""Represents a dependency graph

	Attributes:
		dependencies: Mapping from node to that node's dependents
	"""

	def __init__(self):
		"""Initialize dependencies to empty dictionary"""
		self.dependencies = {}  # type: Dict[str, List[str]]

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
		"""Add an edge such that a depends on b

		If a or b does not exist yet as a node, it will be created.

		Arguments:
			a: The name of the node that depends on b
			b: The name of the node that is depended-upon by a
		"""
		self.dependencies.setdefault(a, []).append(b)
		if b not in self.dependencies:
			self.dependencies[b] = []

	def get_topological_ordering(self):
		# type: () -> List[str]
		"""Get a topological ordering of the nodes

		Returns:
			List of dependency names such that the dependencies can be
			loaded in the order in which they appear in the list without
			violating dependency relationships.

		Raises:
			InvalidDependencyGraphError: If the graph contains a cycle
		"""
		explored = {
			name: ExplorationStatus.UNEXPLORED for name in self.dependencies
		}
		reverse_ordering = []  # type: List[str]
		for node in self.dependencies:
			if explored[node] != ExplorationStatus.EXPLORED:
				self._topo_sort_dfs(node, explored, reverse_ordering)
		return reverse_ordering

	def _topo_sort_dfs(self, node, explored, reverse_ordering):
		# type: (str, Dict[str, bool], List[str]) -> None
		explored[node] = ExplorationStatus.EXPLORING
		for dependency in self.dependencies[node]:
			if explored[dependency] == ExplorationStatus.UNEXPLORED:
				self._topo_sort_dfs(dependency, explored, reverse_ordering)
			elif explored[dependency] == ExplorationStatus.EXPLORING:
				# We have reached a node that we have already
				# visited, but that we have not added to the
				# ordering. Thus a call to dependency is still on
				# the stack. This means there exists a path from
				# dependency to node, so there is a cycle.
				raise InvalidDependencyGraphError(
					"Dependency graphs must not contain cycles")
		reverse_ordering.append(node)
		explored[node] = ExplorationStatus.EXPLORED
