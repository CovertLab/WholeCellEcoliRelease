#!/usr/bin/env python
"""
Constructs a causal network of simulation components along with dynamics data
from a multi-generation simulation.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/26/2018
"""
from __future__ import division

import argparse
import os
import cPickle

import numpy as np
from itertools import izip

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from wholecell.utils.filepath import makedirs

NODE_LIST_HEADER = "ID,class,category,name,synonyms,constants\n"
EDGE_LIST_HEADER = "src_node_id,dst_node_id,stoichiometry,process\n"
DYNAMICS_HEADER = "node,type,units,dynamics\n"

CHECK_SANITY = True
N_GENS = 2

class Node:
	"""
	Class definition for a node in the causal network.
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
		self.dynamics = None
		self.dynamics_units = None

	def get_attributes(self, node_id, name, synonyms="", constants=""):
		"""
		Sets the remaining attribute variables of the node. Argument can be
		in the form of a single dictionary with names of each argument names as
		keys.
		"""
		self.node_id = node_id
		self.name = name
		self.synonyms = synonyms
		self.constants = constants

	def get_dynamics(self, dynamics, dynamics_units):
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
		node_row = "%s,%s,%s,%s,%s,%s\n" % (
			self.node_id, self.node_class, self.node_category,
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

			# Format single string with dynamic attributes separated by commas
			dynamics_row = "%s,%s,%s,%s\n" % (
				self.node_id, name, unit, dynamics
				)

			# Write line to dynamics file
			dynamics_file.write(dynamics_row)


class Edge:
	"""
	Class definition for an edge in the causal network.
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

	def get_attributes(self, src_id, dst_id, stoichiometry=""):
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
		edge_row = "%s,%s,%s,%s\n" % (
			self.src_id, self.dst_id, self.stoichiometry, self.process,
			)

		# Write line to edgelist file
		edgelist_file.write(edge_row)


def add_gene_nodes(simData, simOutDirs, node_list):
	"""
	Add gene nodes with dynamics data to the node list. - Heejo
	"""
	pass

def add_transcript_nodes(simData, simOutDirs, node_list):
	"""
	Add transcript nodes with dynamics data to the node list. - Heejo
	"""
	pass

def add_protein_and_complex_nodes(simData, simOutDirs, node_list):
	"""
	Add protein and complex nodes with dynamics data to the node list. - Eran
	"""
	pass

def add_metabolite_nodes(simData, simOutDirs, node_list):
	"""
	Add metabolite nodes with dynamics data to the node list. - Gwanggyu
	"""
	# TODO: complete this function as a working example
	pass


def add_replication_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add replication nodes with dynamics data to the node list, and add edges
	connected to the replication nodes to the edge list. - Heejo
	"""
	pass

def add_transcription_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add transcription nodes with dynamics data to the node list, and add edges
	connected to the transcription nodes to the edge list. - Heejo
	"""
	pass

def add_translation_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add translation nodes with dynamics data to the node list, and add edges
	connected to the translation nodes to the edge list. - Heejo
	"""
	pass

def add_complexation_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add complexation nodes with dynamics data to the node list, and add edges
	connected to the complexation nodes to the edge list. - Eran
	"""
	pass

def add_metabolism_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add metabolism nodes with dynamics data to the node list, and add edges
	connected to the metabolism nodes to the edge list. - Gwanggyu
	"""
	# TODO: complete this function as a working example
	pass

def add_equilibrium_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add equilibrium nodes with dynamics data to the node list, and add edges
	connected to the equilibrium nodes to the edge list. - Gwanggyu
	"""
	pass

def add_regulation_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add regulation nodes with dynamics data to the node list, and add edges
	connected to the regulation nodes to the edge list. - Gwanggyu
	"""
	pass


def add_time_data(simOutDirs, dynamics_file):
	"""
	Add time data to the dynamics file. - Gwanggyu
	"""
	time = []

	# Loop through all generations
	for simOutDir in simOutDirs:
		# Extract time data from each generation
		t = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		time += list(t)

	time_row = "%s,%s,%s,%s\n" % (
		"time", "time", "s", time
	)

	# Write line to dynamics file
	dynamics_file.write(time_row)


def add_global_dynamics(simData, simOutDirs, dynamics_file):
	"""
	Add global dynamics data to the dynamics file.
	Currently, cell mass and cell volume are the only global dynamics data that
	are being written. - Gwanggyu
	"""
	# Initialize global dynamics lists
	global_dynamics = dict()
	global_dynamics['cell_mass'] = []
	global_dynamics['cell_volume'] = []

	global_dynamics_units = dict()
	global_dynamics_units['cell_mass'] = 'fg'
	global_dynamics_units['cell_volume'] = 'L'

	# Loop through all generations
	for simOutDir in simOutDirs:
		# Extract dynamics data from each generation
		cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
		cell_volume = ((1.0/simData.constants.cellDensity)*(units.fg*cell_mass)).asNumber(units.L)

		# Append to existing list
		global_dynamics['cell_mass'] += list(cell_mass)
		global_dynamics['cell_volume'] += list(cell_volume)

	# Iterate through all dynamics variables associated with the node
	for name, dynamics in global_dynamics.items():
		unit = global_dynamics_units.get(name, "")

		# Format single string with dynamic attributes separated by commas
		dynamics_row = "%s,%s,%s,%s\n" % (
			"global", name, unit, dynamics
		)

		# Write line to dynamics file
		dynamics_file.write(dynamics_row)


def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile=None, metadata=None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
	assert ap.n_generation >= N_GENS

	simData = cPickle.load(open(simDataFile))

	# Get first cell from each generation
	first_cell_lineage = []

	# For all generation indexes subject to analysis, get first cell
	for gen_idx in range(N_GENS):
		first_cell_lineage.append(ap.get_cells(generation=[gen_idx])[0])

	simOutDirs = []

	# Go through first cells in each generation
	for gen, simDir in enumerate(first_cell_lineage):
		simOutDir = os.path.join(simDir, "simOut")
		simOutDirs.append(simOutDir)

	# Initialize node list and edge list
	node_list = []
	edge_list = []

	# Add state nodes to the node list
	add_gene_nodes(simData, simOutDirs, node_list)  # Heejo
	add_transcript_nodes(simData, simOutDirs, node_list)  # Heejo
	add_protein_and_complex_nodes(simData, simOutDirs, node_list)  # Eran
	add_metabolite_nodes(simData, simOutDirs, node_list)  # Gwanggyu

	# Add process nodes and associated edges to the node list and edge list, respectively
	add_replication_nodes_and_edges(simData, simOutDirs, node_list, edge_list)  # Heejo
	add_transcription_nodes_and_edges(simData, simOutDirs, node_list, edge_list)  # Heejo
	add_translation_nodes_and_edges(simData, simOutDirs, node_list, edge_list)  # Heejo
	add_complexation_nodes_and_edges(simData, simOutDirs, node_list, edge_list)  # Eran
	add_metabolism_nodes_and_edges(simData, simOutDirs, node_list, edge_list)  # Gwanggyu
	add_equilibrium_nodes_and_edges(simData, simOutDirs, node_list, edge_list)  # Gwanggyu
	add_regulation_nodes_and_edges(simData, simOutDirs, node_list, edge_list)  # Gwanggyu

	# TODO: Check for network sanity (optional)
	if CHECK_SANITY:
		pass

	# Open node/edge list files and dynamics file
	nodelist_file = open(os.path.join(plotOutDir, plotOutFileName + "_nodelist.csv"), 'w')
	edgelist_file = open(os.path.join(plotOutDir, plotOutFileName + "_edgelist.csv"), 'w')
	dynamics_file = open(os.path.join(plotOutDir, plotOutFileName + "_dynamics.csv"), 'w')

	# Write header rows to each of the files
	nodelist_file.write(NODE_LIST_HEADER)
	edgelist_file.write(EDGE_LIST_HEADER)
	dynamics_file.write(DYNAMICS_HEADER)

	# Add time and global dynamics data to dynamics file
	add_time_data(simOutDirs, dynamics_file)
	add_global_dynamics(simData, simOutDirs, dynamics_file)

	# Write node, edge list and dynamics data csv files
	for node in node_list:
		node.write_nodelist(nodelist_file)
		node.write_dynamics(dynamics_file)

	for edge in edge_list:
		edge.write_edgelist(edgelist_file)


if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
		wholecell.utils.constants.SERIALIZED_KB_DIR,
		wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
	)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help="Directory containing simulation output", type=str)
	parser.add_argument("plotOutDir", help="Directory containing plot output (will get created if necessary)", type=str)
	parser.add_argument("plotOutFileName", help="File name to produce", type=str)
	parser.add_argument("--simDataFile", help="KB file name", type=str, default=defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
