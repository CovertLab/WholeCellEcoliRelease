#!/usr/bin/env python
"""
Classes for the Nodes and Edges of a causality network.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/28/2018
"""
from __future__ import absolute_import
from __future__ import division

NODELIST_FILENAME = "causality_network_node_list.tsv"
EDGELIST_FILENAME = "causality_network_edge_list.tsv"

DYNAMICS_PRECISION = 6
PROBABILITY_PRECISION = 4
TIME_PRECISION = 2

class Node:
	"""
	Class definition for a node in the causality network.
	"""

	def __init__(self, node_class, node_type):
		"""
		Initializes instance variables. Node class and type must be given as
		arguments.

		Variables:
			node_class: Class of node, string, either "State" or "Process"
			node_type: Type of node, string, e.g. "Gene", "Metabolism"
			node_id: Unique ID of node, string, e.g. "EG11274", "CPLX-125[c]"
			name: Generic name of node, string, e.g. "trpL", "pyruvate"
			synonyms: List of synonyms of node, list of strings
				e.g. ["anth", "tryD", tryp-4"]
			constants: Dictionary with constant names as keys and constants as
				values, dictionary, e.g. {"reversibility": 0, "Km": 1e-6}
			dynamics: Dictionary with dynamics data type as keys and list of
				time-series data as values, dictionary,
				e.g. {"counts": [8151, 8525, ...],
					  "concentration": [1.151e-7, 1.155e-7, ...]}
			dynamics_units: Dictionary with dynamics data type as keys and its
				units as values (must share same keys with dynamics),
				dictionary, e.g. {"counts": "N", "concentration": "mol/L"}
		"""
		self.node_class = node_class
		self.node_type = node_type
		self.node_id = None
		self.name = None
		self.synonyms = None
		self.constants = None
		self.dynamics = {}
		self.dynamics_units = {}

	def get_node_id(self):
		"""
		Return ID of node.
		"""
		return self.node_id

	def read_attributes(self, node_id, name, synonyms="", constants=""):
		"""
		Sets the remaining attribute variables of the node. Argument can be
		in the form of a single dictionary with names of each argument names as
		keys.
		"""
		self.node_id = node_id
		self.name = name
		self.synonyms = synonyms
		self.constants = constants

	def read_dynamics(self, dynamics, dynamics_units):
		"""
		Sets the dynamics variable of the node.
		"""
		self.dynamics = dynamics
		self.dynamics_units = dynamics_units

	def write_nodelist(self, nodelist_file):
		"""
		Writes a single row specifying the given node to the nodelist file.
		"""
		# Format single string with attributes of the node separated by commas
		node_row = "%s\t%s\t%s\t%s\t%s\t%s\n" % (
			self.node_id, self.node_class, self.node_type,
			self.name, self.synonyms, self.constants,
			)

		# Write line to nodelist file
		nodelist_file.write(node_row)

	def write_dynamics(self, dynamics_file):
		"""
		Writes a single row of dynamics data for each dynamics variable
		associated with the node.
		"""
		# Iterate through all dynamics variables associated with the node
		for name, dynamics in self.dynamics.items():
			unit = self.dynamics_units.get(name, "")

			# Format dynamics string depending on data type
			if unit == "N":
				dynamics_string = self._format_dynamics_string(dynamics, "int")
			elif unit == "prob":
				dynamics_string = self._format_dynamics_string(dynamics, "prob")
			else:
				dynamics_string = self._format_dynamics_string(dynamics, "float")

			# Format single string with dynamic attributes separated by commas
			dynamics_row = "%s\t%s\t%s\t%s\n" % (
				self.node_id, name, unit, dynamics_string
				)

			# Write line to dynamics file
			dynamics_file.write(dynamics_row)

	def _format_dynamics_string(self, dynamics, datatype):
		"""
		Formats the string of dynamics data that is printed out to the dynamics
		file. If datatype is "int", print all numbers as full decimal integers.
		If datatype is "float", print all numbers in the general format with the
		precision set by DYNAMICS_PRECISION. If datatype is "time", print all
		numbers in the floating point format with the precision set by
		TIME_PRECISION.
		"""
		if datatype == "int":
			dynamics_string = ", ".join(
				"{0:d}".format(val) for val in dynamics)

		elif datatype == "float":
			dynamics_string = ", ".join(
				"{0:.{1}g}".format(val, DYNAMICS_PRECISION) for val in
				dynamics)

		elif datatype == "prob":
			dynamics_string = ", ".join(
				"{0:.{1}f}".format(val, PROBABILITY_PRECISION) for val in
				dynamics)

		elif datatype == "time":
			dynamics_string = ", ".join(
				"{0:.{1}f}".format(val, TIME_PRECISION) for val in dynamics)

		else:
			dynamics_string = dynamics

		return dynamics_string


class Edge:
	"""
	Class definition for an edge in the causality network.
	"""

	def __init__(self, process):
		"""
		Initializes instance variables. Edge type must be given as arguments.

		Variables:
			edge_type: Type of edge (type of the process node the edge is
				attached to), string, e.g. "Complexation", "Metabolism"
			src_id: ID of the source node, string, e.g. "RXN0-2382"
			dst_id: ID of the destination node, string, e.g. "WATER[c]"
			stoichiometry: (Only for metabolism edges) Stoichiometric
				coefficient of reaction-metabolite pair, integer, e.g. 1
		"""
		self.process = process
		self.src_id = None
		self.dst_id = None
		self.stoichiometry = None

	def get_src_id(self):
		"""
		Return ID of source node.
		"""
		return self.src_id

	def get_dst_id(self):
		"""
		Return ID of destination node.
		"""
		return self.dst_id

	def get_process(self):
		"""
		Return process associated with the edge.
		"""
		return self.process

	def read_attributes(self, src_id, dst_id, stoichiometry=""):
		"""
		Sets the remaining attribute variables of the node. Argument can be
		in the form of a single dictionary with names of each argument names as
		keys.
		"""
		self.src_id = src_id
		self.dst_id = dst_id
		self.stoichiometry = stoichiometry

	def write_edgelist(self, edgelist_file):
		"""
		Writes a single row specifying the given edge to the edgelist file.
		"""
		# Format single string with attributes of the edge separated by commas
		edge_row = "%s\t%s\t%s\t%s\n" % (
			self.src_id, self.dst_id, self.stoichiometry, self.process,
			)

		# Write line to edgelist file
		edgelist_file.write(edge_row)
