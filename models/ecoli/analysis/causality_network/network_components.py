"""
Classes for the Nodes and Edges of a causality network.
"""
from __future__ import absolute_import, division, print_function

from typing import Optional, Union
import six


# Filenames
NODELIST_FILENAME = "causality_network_node_list.tsv"
EDGELIST_FILENAME = "causality_network_edge_list.tsv"
DYNAMICS_FILENAME = "causality_network_dynamics.tsv"
NODELIST_JSON = 'nodes.json'
EDGELIST_JSON = 'edges.json'

# Headers
NODE_LIST_HEADER = "\t".join(
	["ID", "class", "type", "name", "synonyms", "constants", "url"]
	)
EDGE_LIST_HEADER = "\t".join(
	["src_node_id", "dst_node_id", "stoichiometry", "process"]
	)
DYNAMICS_HEADER = "\t".join(
	["ID", "type", "units", "dynamics"]
	)

# Special strings used as units to designate type of dynamics
COUNT_UNITS = "N"
PROB_UNITS = "prob"

# Precision settings for numbers in the dynamics file
DYNAMICS_PRECISION = 6
PROBABILITY_PRECISION = 4
TIME_PRECISION = 2

class Node(object):
	"""
	Class definition for a node in the causality network.
	"""

	def __init__(self):
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
			url: URL to EcoCyc page, string, eg. "https://ecocyc.org/ECOLI/
				substring-search?type=NIL&object=EG11028&quickSearch=Quick+
				Search"
		"""
		self.node_class = None
		self.node_type = None
		self.node_id = None
		self.name = None
		self.synonyms = None
		self.constants = None
		self.dynamics = {}
		self.dynamics_units = {}
		self.url = None
		self.location = None

	def get_node_id(self):
		"""
		Return ID of node.
		"""
		return self.node_id

	def read_attributes(self, node_class, node_type, node_id, name="",
						synonyms="", constants="", url="", location=""):
		"""
		Sets the attribute variables of the node. Argument can be in the form
		of a single dictionary with names of each argument names as keys.
		"""
		self.node_class = node_class
		self.node_type = node_type
		self.node_id = node_id
		self.name = name
		self.synonyms = synonyms
		self.constants = constants
		self.url = url
		self.location = location

	def read_attributes_from_tsv(self, tsv_line):
		"""
		Reads attributes (node type and node id) from a tab-delimited line in
		the node_list.tsv file.
		"""
		split_tsv_line = tsv_line[:-1].split('\t')

		self.node_type = split_tsv_line[2]
		self.node_id = split_tsv_line[0]

		return (self.node_id, self.node_type)


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
		node_row = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
			self.node_id, self.node_class, self.node_type,
			self.name, self.synonyms, self.constants, self.url,
			)

		# Write line to nodelist file
		nodelist_file.write(node_row + "\n")

	def write_dynamics(self, dynamics_file):
		"""
		Writes a single row of dynamics data for each dynamics variable
		associated with the node.
		"""
		# Iterate through all dynamics variables associated with the node
		for dynamics_name, dynamics_data in six.viewitems(self.dynamics):
			unit = self.dynamics_units.get(dynamics_name, "")

			# Format dynamics string depending on data type
			if unit == COUNT_UNITS:
				dynamics_string = self._format_dynamics_string(dynamics_data, "int")
			elif unit == PROB_UNITS:
				dynamics_string = self._format_dynamics_string(dynamics_data, "prob")
			elif unit == "s":
				dynamics_string = self._format_dynamics_string(dynamics_data, "time")
			else:
				dynamics_string = self._format_dynamics_string(dynamics_data, "float")

			# Format single string with dynamic attributes separated by commas
			dynamics_row = "%s\t%s\t%s\t%s" % (
				self.node_id, dynamics_name, unit, dynamics_string
				)

			# Write line to dynamics file
			dynamics_file.write(dynamics_row + "\n")

	def dynamics_dict(self):
		all_dynamics = []
		for name, data in six.viewitems(self.dynamics):
			unit = self.dynamics_units.get(name, "")
			dynamics = {
				'units': unit,
				'type': name,
				'id': self.node_id,
				'dynamics': data.tolist()}
			all_dynamics.append(dynamics)
		return all_dynamics

	def to_dict(self):
		synonyms = []
		if isinstance(self.synonyms, list):
			synonyms = self.synonyms

			# Some of the synonyms are strings of list-like entities and some are actual lists
			# --------------------------------------------------------------------------------
			# try:
			# 	synonyms = ast.literal_eval(self.synonyms or '[]')
			# except:
			# 	print('parsing synonyms failed for {} {}'.format(type(self.synonyms), self.synonyms))

		return {
			'ID': self.node_id,
			'type': self.node_type,
			'name': self.name,
			'class': self.node_class,
			'synonyms': synonyms,
			'constants': self.constants,
			'url': self.url,
			'location': self.location}


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


class Edge(object):
	"""
	Class definition for an edge in the causality network.
	"""

	def __init__(self, process):
		# type: (str) -> None
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
		self.process = process  # type: str
		self.src_id = None  # type: Optional[str]
		self.dst_id = None  # type: Optional[str]
		self.stoichiometry = None  # type: Union[None, str, int]

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
		# type: (str, str, Union[None, str, int]) -> None
		# TODO(jerry): A narrower type for stoichiometry?
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
		edge_row = "%s\t%s\t%s\t%s" % (
			self.src_id, self.dst_id, self.stoichiometry, self.process,
			)

		# Write line to edgelist file
		edgelist_file.write(edge_row + "\n")
