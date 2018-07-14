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
DYNAMICS_PRECISION = 6
TIME_PRECISION = 2

# Proteins that are reactants and products of a metabolic reaction
PROTEINS_IN_METABOLISM = ["EG50003-MONOMER[c]", "PHOB-MONOMER[c]", "PTSI-MONOMER[c]", "PTSH-MONOMER[c]"]

# Equilibrium complexes that are formed from deleted equilibrium reactions, but
# are reactants in a complexation reaction
EQUILIBRIUM_COMPLEXES_IN_COMPLEXATION = ["CPLX0-7620[c]", "CPLX0-7701[c]", "CPLX0-7677[c]", "MONOMER0-1781[c]", "CPLX0-7702[c]"]

# Metabolites that are used as ligands in equilibrium, but do not participate
# in any metabolic reactions
METABOLITES_ONLY_IN_EQUILIBRIUM = ["4FE-4S[c]", "NITRATE[p]"]

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
		node_row = "%s,%s,%s,%s,%s,%s\n" % (
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
				dynamics_string = format_dynamics_string(dynamics, "int")
			else:
				dynamics_string = format_dynamics_string(dynamics, "float")

			# Format single string with dynamic attributes separated by commas
			dynamics_row = "%s,%s,%s,%s\n" % (
				self.node_id, name, unit, dynamics_string
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
		edge_row = "%s,%s,%s,%s\n" % (
			self.src_id, self.dst_id, self.stoichiometry, self.process,
			)

		# Write line to edgelist file
		edgelist_file.write(edge_row)


def add_gene_nodes(simData, simOutDirs, node_list):
	"""
	Add gene nodes with dynamics data to the node list. - Heejo
	"""
	# This is currently being done by the add_replication_nodes_and_edges()
	# function.
	pass

def add_transcript_nodes(simData, simOutDirs, node_list):
	"""
	Add transcript nodes with dynamics data to the node list. - Heejo
	"""
	# This is currently being done by the add_transcription_nodes_and_edges()
	# function.
	pass

def add_protein_and_complex_nodes(simData, simOutDirs, node_list):
	"""
	Add protein and complex nodes with dynamics data to the node list. - Eran
	"""
	# This is currently being done by the add_complexation_nodes_and_edges()
	# function.
	pass

def add_metabolite_nodes(simData, simOutDirs, node_list):
	"""
	Add metabolite nodes with dynamics data to the node list. - Gwanggyu
	"""
	# This is currently being done by the add_metabolism_nodes_and_edges()
	# function.
	pass

def add_replication_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add replication nodes with dynamics data to the node list, and add edges
	connected to the replication nodes to the edge list. - Heejo
	"""
	dntp_ids = simData.moleculeGroups.dNtpIds
	ppi_id = "PPI[c]"
	dnap_ids = ['CPLX0-2361[c]', 'CPLX0-3761[c]', 'CPLX0-3925[c]', 'CPLX0-7910[c]']

	# Loop through all genes
	geneIds = [data[0] for data in simData.process.replication.geneData]

	for geneId in geneIds:
		# Initialize a single gene node
		gene_node = Node("State", "Gene")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': geneId, 'name': geneId}
		gene_node.read_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append gene node to node_list
		node_list.append(gene_node)

		# Initialize a single replication node for each gene
		replication_node = Node("Process", "Replication")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		replication_node_id = "%s_REPLICATION" % geneId
		attr = {'node_id': replication_node_id, 'name': replication_node_id}
		replication_node.read_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append replication node to node_list
		node_list.append(replication_node)

		# Add edges between gene and replication nodes
		gene_to_replication_edge = Edge("Replication")
		attr = {'src_id': geneId, 'dst_id': replication_node_id}
		gene_to_replication_edge.read_attributes(**attr)
		edge_list.append(gene_to_replication_edge)

		replication_to_gene_edge = Edge("Replication")
		attr = {'src_id': replication_node_id, 'dst_id': geneId}
		replication_to_gene_edge.read_attributes(**attr)
		edge_list.append(replication_to_gene_edge)

		# Add edges from dNTPs to replication nodes
		for dntp_id in dntp_ids:
			dNTP_to_replication_edge = Edge("Replication")
			attr = {'src_id': dntp_id, 'dst_id': replication_node_id}
			dNTP_to_replication_edge.read_attributes(**attr)
			edge_list.append(dNTP_to_replication_edge)

		# Add edge from replication to Ppi
		replication_to_ppi_edge = Edge("Replication")
		attr = {'src_id': replication_node_id, 'dst_id': ppi_id}
		replication_to_ppi_edge.read_attributes(**attr)
		edge_list.append(replication_to_ppi_edge)

		# Add edges from DNA polymerases to replication
		for dnap_id in dnap_ids:
			pol_to_replication_edge = Edge("Replication")
			attr = {'src_id': dnap_id, 'dst_id': replication_node_id}
			pol_to_replication_edge.read_attributes(**attr)
			edge_list.append(pol_to_replication_edge)


def add_transcription_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add transcription nodes with dynamics data to the node list, and add edges
	connected to the transcription nodes to the edge list. - Heejo
	"""
	ntp_ids = simData.moleculeGroups.ntpIds
	ppi_id = "PPI[c]"
	rnap_id = "APORNAP-CPLX[c]"

	# Loop through all genes
	for geneId, rnaId, _ in simData.process.replication.geneData:
		# Initialize a single transcript node
		rna_node = Node("State", "RNA")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		rna_node_id = "%s[c]" % rnaId
		attr = {'node_id': rna_node_id, 'name': rna_node_id}
		rna_node.read_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append transcript node to node_list
		node_list.append(rna_node)

		# Initialize a single transcription node for each gene
		transcription_node = Node("Process", "Transcription")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		transcription_node_id = "%s_TRANSCRIPTION" % geneId
		attr = {'node_id': transcription_node_id, 'name': transcription_node_id}
		transcription_node.read_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append transcription node to node_list
		node_list.append(transcription_node)

		# Add edge from gene to transcription node
		gene_to_transcription_edge = Edge("Transcription")
		attr = {'src_id': geneId, 'dst_id': transcription_node_id}
		gene_to_transcription_edge.read_attributes(**attr)
		edge_list.append(gene_to_transcription_edge)

		# Add edge from transcription to transcript node
		transcription_to_rna_edge = Edge("Transcription")
		attr = {'src_id': transcription_node_id, 'dst_id': rna_node_id}
		transcription_to_rna_edge.read_attributes(**attr)
		edge_list.append(transcription_to_rna_edge)

		# Add edges from NTPs to transcription nodes
		for ntp_id in ntp_ids:
			ntp_to_transcription_edge = Edge("Transcription")
			attr = {'src_id': ntp_id, 'dst_id': transcription_node_id}
			ntp_to_transcription_edge.read_attributes(**attr)
			edge_list.append(ntp_to_transcription_edge)

		# Add edge from transcription to Ppi
		transcription_to_ppi_edge = Edge("Transcription")
		attr = {'src_id': transcription_node_id, 'dst_id': ppi_id}
		transcription_to_ppi_edge.read_attributes(**attr)
		edge_list.append(transcription_to_ppi_edge)

		# Add edges from RNA polymerases to transcription
		pol_to_transcription_edge = Edge("Transcription")
		attr = {'src_id': rnap_id, 'dst_id': transcription_node_id}
		pol_to_transcription_edge.read_attributes(**attr)
		edge_list.append(pol_to_transcription_edge)


def add_translation_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add translation nodes with dynamics data to the node list, and add edges
	connected to the translation nodes to the edge list. - Heejo
	"""
	# Create nodes for amino acids
	aa_ids = simData.moleculeGroups.aaIDs
	gtp_id = "GTP[c]"
	gdp_id = "GDP[c]"
	water_id = "WATER[c]"
	ppi_id = "PPI[c]"

	ribosome_subunit_ids = [simData.moleculeGroups.s30_fullComplex[0], simData.moleculeGroups.s50_fullComplex[0]]

	# Loop through all translatable genes
	for data in simData.process.translation.monomerData:
		monomerId = data[0]
		rnaId = data[1]
		geneId = rnaId.split("_RNA[c]")[0]

		# Initialize a single protein node
		protein_node = Node("State", "Protein")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': monomerId, 'name': monomerId}
		protein_node.read_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append protein node to node_list
		node_list.append(protein_node)

		# Initialize a single translation node for each transcript
		translation_node = Node("Process", "Translation")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		translation_node_id = "%s_TRANSLATION" % geneId
		attr = {'node_id': translation_node_id, 'name': translation_node_id}
		translation_node.read_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append translation node to node_list
		node_list.append(translation_node)

		# Add edge from transcript to translation node
		rna_to_translation_edge = Edge("Translation")
		attr = {'src_id': rnaId, 'dst_id': translation_node_id}
		rna_to_translation_edge.read_attributes(**attr)
		edge_list.append(rna_to_translation_edge)

		# Add edge from translation to monomer node
		translation_to_protein_edge = Edge("Translation")
		attr = {'src_id': translation_node_id, 'dst_id': monomerId}
		translation_to_protein_edge.read_attributes(**attr)
		edge_list.append(translation_to_protein_edge)

		# Add edges from amino acids to translation node
		for aa_id in aa_ids:
			aa_to_translation_edge = Edge("Translation")
			attr = {'src_id': aa_id, 'dst_id': translation_node_id}
			aa_to_translation_edge.read_attributes(**attr)
			edge_list.append(aa_to_translation_edge)

		# Add edges from other reactants to translation node
		for reactant_id in [gtp_id, water_id]:
			reactant_to_translation_edge = Edge("Translation")
			attr = {'src_id': reactant_id, 'dst_id': translation_node_id}
			reactant_to_translation_edge.read_attributes(**attr)
			edge_list.append(reactant_to_translation_edge)

		# Add edges from translation to other product nodes
		for product_id in [gdp_id, ppi_id, water_id]:
			translation_to_product_edge = Edge("Translation")
			attr = {'src_id': translation_node_id, 'dst_id': product_id}
			translation_to_product_edge.read_attributes(**attr)
			edge_list.append(translation_to_product_edge)

		# Add edges from ribosome subunits to translation node
		for subunit_id in ribosome_subunit_ids:
			subunit_to_translation_edge = Edge("Translation")
			attr = {'src_id': subunit_id, 'dst_id': translation_node_id}
			subunit_to_translation_edge.read_attributes(**attr)
			edge_list.append(subunit_to_translation_edge)


def add_complexation_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add complexation nodes with dynamics data to the node list, and add edges
	connected to the complexation nodes to the edge list. - Eran
	"""
	simOutDir = simOutDirs[0]

	# TODO (Eran) raw_data is here used to get complexation reaction IDs and stoichiometry. This can be saved to sim_data and then retrieved here.
	from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
	raw_data = KnowledgeBaseEcoli()

	# get reaction IDs from raw data
	reactionIDs = [dict['id'] for dict in raw_data.complexationReactions]

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get dynamics data from all simOutDirs (# rxns/ts for complexation, counts
	# for proteins and complexes)
	reactions_array = np.empty((0, len(reactionIDs)))
	counts_array = np.empty((0, len(moleculeIDs)))

	for simOutDir in simOutDirs:
		#TODO (ERAN) save complex reaction rate (# rxns/ts) in reaction_array. Save this in listener, or compute it here.
		# fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		# reactionFluxes = fbaResults.readColumn('reactionFluxes')
		# reactions_array = np.concatenate((reactions_array, reactionFluxes))

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

	# Get complexation stoichiometry from simData
	complexStoich = {}
	for reaction in raw_data.complexationReactions:
		stoich = {}
		for molecule in reaction['stoichiometry']:
			molecule_name = '%s[%s]' % (molecule['molecule'], molecule['location'])
			stoich[molecule_name] = molecule['coeff']
		complexStoich[reaction['id']] = stoich

	# List of all complex IDs
	complex_ids = simData.process.complexation.ids_complexes + EQUILIBRIUM_COMPLEXES_IN_COMPLEXATION

	# Loop through all complexation reactions
	for idx, reaction in enumerate(reactionIDs):
		# Initialize a single complexation node for each complexation reaction
		complexation_node = Node("Process", "Complexation")

		# Add attributes to the node
		attr = {'node_id': reaction, 'name': reaction}
		complexation_node.read_attributes(**attr)

		# # TODO (ERAN) Add dynamics data (# rxns/ts) to the node.
		# dynamics = {'flux': list(flux_array[:, idx])}
		# dynamics_units = {'flux': 'mmol/gCDW/h'}
		# complexation_node.read_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(complexation_node)

		# Get reaction stoichiometry from complexStoich
		stoich_dict = complexStoich[reaction]

		# Loop through all proteins participating in the reaction
		for protein, stoich in stoich_dict.items():

			# Initialize complex edge
			complex_edge = Edge("Complexation")

			# Add attributes to the complex edge
			# Note: the direction of the edge is determined by the sign of the
			# stoichiometric coefficient.
			if stoich > 0:
				attr = {'src_id': reaction,
					'dst_id': protein,
					'stoichiometry': stoich
					}
			else:
				attr = {'src_id': protein,
					'dst_id': reaction,
					'stoichiometry': stoich
					}
			complex_edge.read_attributes(**attr)

			# Append edge to edge_list
			edge_list.append(complex_edge)

	for complex in complex_ids:
		# Initialize a single complex node for each complex
		complex_node = Node("State", "Complex")

		# Add attributes to the node
		# TODO: Get molecular mass using getMass().
		# TODO: Get correct protein name and synonyms from EcoCyc
		attr = {'node_id': complex,
			'name': complex,
			'constants': {'mass': 0}
			}
		complex_node.read_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the complex in the counts array
		try:
			complex_idx = moleculeIDs.index(complex)
		except ValueError:  # complex ID not found in moleculeIDs
			complex_idx = -1

		if complex_idx != -1:
			dynamics = {'counts': list(counts_array[:, complex_idx].astype(np.int))}
			dynamics_units = {'counts': 'N'}
			complex_node.read_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(complex_node)


def add_metabolism_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add metabolism nodes with dynamics data to the node list, and add edges
	connected to the metabolism nodes to the edge list. - Gwanggyu
	Note: forward and reverse reactions are represented as separate nodes.
	State nodes for metabolites are also added here.
	"""
	# Get reaction list from first simOut directory
	simOutDir = simOutDirs[0]
	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = fbaResults.readAttribute("reactionIDs")

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get dynamics data from all simOutDirs (flux for reactions, counts for
	# metabolites)
	flux_array = np.empty((0, len(reactionIDs)))
	counts_array = np.empty((0, len(moleculeIDs)), dtype=np.int)

	for simOutDir in simOutDirs:
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionFluxes = fbaResults.readColumn('reactionFluxes')
		flux_array = np.concatenate((flux_array, reactionFluxes))

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

	# Get all reaction stoichiometry from simData
	reactionStoich = simData.process.metabolism.reactionStoich

	# Get reaction to catalyst dict from simData
	reactionCatalysts = simData.process.metabolism.reactionCatalysts

	# Initialize list of metabolite IDs
	metabolite_ids = []

	# Loop through all reactions
	for idx, reaction in enumerate(reactionIDs):
		# Initialize a single metabolism node for each reaction
		metabolism_node = Node("Process", "Metabolism")

		# Add attributes to the node
		# TODO: Get correct reaction name and synonyms from EcoCyc
		attr = {'node_id': reaction, 'name': reaction}
		metabolism_node.read_attributes(**attr)

		# Add dynamics data (flux) to the node. The flux array shares the same
		# column index with the reactionIDs list
		dynamics = {'flux': list(flux_array[:, idx])}
		dynamics_units = {'flux': 'mmol/gCDW/h'}
		metabolism_node.read_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(metabolism_node)

		# Get reaction stoichiometry from reactionStoich
		stoich_dict = reactionStoich[reaction]

		# Get list of proteins that catalyze this reaction
		catalyst_list = reactionCatalysts.get(reaction, [])

		# Add an edge from each catalyst to the metabolism node
		for catalyst in catalyst_list:
			metabolism_edge = Edge("Metabolism")
			attr = {'src_id': catalyst,
				'dst_id': reaction
			}

			metabolism_edge.read_attributes(**attr)
			edge_list.append(metabolism_edge)

		# Loop through all metabolites participating in the reaction
		for metabolite, stoich in stoich_dict.items():
			# Add metabolites that were not encountered
			if metabolite not in metabolite_ids:
				metabolite_ids.append(metabolite)

			# Initialize Metabolism edge
			metabolism_edge = Edge("Metabolism")

			# Add attributes to the Metabolism edge
			# Note: the direction of the edge is determined by the sign of the
			# stoichiometric coefficient.
			if stoich > 0:
				attr = {'src_id': reaction,
					'dst_id': metabolite,
					'stoichiometry': stoich
					}
			else:
				attr = {'src_id': metabolite,
					'dst_id': reaction,
					'stoichiometry': stoich
					}
			metabolism_edge.read_attributes(**attr)

			# Append edge to edge_list
			edge_list.append(metabolism_edge)

	# Loop through all metabolites
	for metabolite in metabolite_ids:
		if metabolite in PROTEINS_IN_METABOLISM:
			continue

		# Initialize a single metabolite node for each metabolite
		metabolite_node = Node("State", "Metabolite")

		# Add attributes to the node
		# TODO: Get molecular mass using getMass(). Some of the metabolites do not have mass data?
		# TODO: Get correct metabolite name and synonyms from EcoCyc
		attr = {'node_id': metabolite,
			'name': metabolite,
			'constants': {'mass': 0}
			}
		metabolite_node.read_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the metabolite in the counts array
		try:
			metabolite_idx = moleculeIDs.index(metabolite)
		except ValueError:  # metabolite ID not found in moleculeIDs
			# Some of the metabolites are not being tracked.
			# TODO: why?
			metabolite_idx = -1

		if metabolite_idx != -1:
			dynamics = {'counts': list(counts_array[:, metabolite_idx].astype(np.int))}
			dynamics_units = {'counts': 'N'}
			metabolite_node.read_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(metabolite_node)


def add_equilibrium_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add equilibrium nodes with dynamics data to the node list, and add edges
	connected to the equilibrium nodes to the edge list. - Gwanggyu
	"""
	# Get equilibrium-specific data from simData
	moleculeIds = simData.process.equilibrium.moleculeNames
	rxnIds = simData.process.equilibrium.rxnIds
	stoichMatrix = simData.process.equilibrium.stoichMatrix()
	ratesFwd = np.array(simData.process.equilibrium.ratesFwd, dtype=np.float32)
	ratesRev = np.array(simData.process.equilibrium.ratesRev, dtype=np.float32)

	# Get bulkMolecule IDs from first simOut directory
	simOutDir = simOutDirs[0]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get dynamics data from all simOutDirs
	counts_array = np.empty((0, len(moleculeIDs)), dtype=np.int)

	for simOutDir in simOutDirs:
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

	# Get IDs of complexes that were already added
	complexation_complex_ids = simData.process.complexation.ids_complexes

	# Get list of complex IDs in equilibrium
	equilibrium_complex_ids = simData.process.equilibrium.ids_complexes

	# TODO: get dynamics data for equilibrium nodes. Will need new listener.

	# Loop through each equilibrium reaction
	for reactionIdx, rxnId in enumerate(rxnIds):

		# Initialize a single equilibrium node for each equilibrium reaction
		equilibrium_node = Node("Process", "Equilibrium")

		# Add attributes to the node
		rxnName = rxnId[:-4] + " equilibrium reaction"
		attr = {'node_id': rxnId,
			'name': rxnName,
			'constants': {'rateFwd': ratesFwd[reactionIdx], 'rateRev': ratesRev[reactionIdx]}
		}
		equilibrium_node.read_attributes(**attr)

		# Append new node to node_list
		node_list.append(equilibrium_node)

		# Extract column corresponding to reaction in the stoichiometric matrix
		stoichMatrixColumn = stoichMatrix[:, reactionIdx]

		# Loop through each element in column
		for moleculeIdx, stoich in enumerate(stoichMatrixColumn):
			moleculeId = moleculeIds[moleculeIdx]

			# If the stoichiometric coefficient is negative, add reactant edge
			# to the equilibrium node
			if stoich < 0:
				equilibrium_edge = Edge("Equilibrium")
				attr = {'src_id': moleculeId,
					'dst_id': rxnId,
					'stoichiometry': stoich
				}

				equilibrium_edge.read_attributes(**attr)
				edge_list.append(equilibrium_edge)

			# If the coefficient is positive, add product edge
			elif stoich > 0:
				equilibrium_edge = Edge("Equilibrium")
				attr = {'src_id': rxnId,
					'dst_id': moleculeId,
					'stoichiometry': stoich
				}

				equilibrium_edge.read_attributes(**attr)
				edge_list.append(equilibrium_edge)

	for complex_id in equilibrium_complex_ids:
		if complex_id in complexation_complex_ids:
			continue

		# Initialize a single complex node for each complex
		complex_node = Node("State", "Complex")

		# Add attributes to the node
		# TODO: Get molecular mass using getMass().
		# TODO: Get correct protein name and synonyms from EcoCyc
		attr = {'node_id': complex_id,
			'name': complex_id,
			'constants': {'mass': 0}
			}
		complex_node.read_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the complex in the counts array
		try:
			complex_idx = moleculeIDs.index(complex_id)
		except ValueError:  # complex ID not found in moleculeIDs
			complex_idx = -1

		if complex_idx != -1:
			dynamics = {'counts': list(counts_array[:, complex_idx].astype(np.int))}
			dynamics_units = {'counts': 'N'}
			complex_node.read_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(complex_node)

	# Loop through all metabolites
	for metabolite in METABOLITES_ONLY_IN_EQUILIBRIUM:
		# Initialize a single metabolite node for each metabolite
		metabolite_node = Node("State", "Metabolite")

		# Add attributes to the node
		# TODO: Get molecular mass using getMass(). Some of the metabolites do not have mass data?
		# TODO: Get correct metabolite name and synonyms from EcoCyc
		attr = {'node_id': metabolite,
			'name': metabolite,
			'constants': {'mass': 0}
			}
		metabolite_node.read_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the metabolite in the counts array
		try:
			metabolite_idx = moleculeIDs.index(metabolite)
		except ValueError:  # metabolite ID not found in moleculeIDs
			# Some of the metabolites are not being tracked.
			# TODO: why?
			metabolite_idx = -1

		if metabolite_idx != -1:
			dynamics = {'counts': list(counts_array[:, metabolite_idx].astype(np.int))}
			dynamics_units = {'counts': 'N'}
			metabolite_node.read_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(metabolite_node)


def add_regulation_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add regulation nodes with dynamics data to the node list, and add edges
	connected to the regulation nodes to the edge list. - Gwanggyu
	"""
	# Get regulation-specific data from simData
	tfToFC = simData.tfToFC

	# Get bulkMolecule IDs from first simOut directory
	simOutDir = simOutDirs[0]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get dynamics data from all simOutDirs
	counts_array = np.empty((0, len(moleculeIDs)), dtype=np.int)

	# TF-DNA bound counts are stored in bulkMolecules counts
	for simOutDir in simOutDirs:
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

	# Get IDs of genes and RNAs
	gene_ids = []
	rna_ids = []

	for gene_id, rna_id, _ in simData.process.replication.geneData:
		gene_ids.append(gene_id)
		rna_ids.append(rna_id)

	# Loop through all TFs
	for tf, transcriptIDdict in tfToFC.items():
		# Add localization ID to the TF ID
		tfID = tf + "[c]"

		# Loop through all transcripts the TF regulates
		for transcriptID in transcriptIDdict.keys():

			# Get ID of TF-gene pair
			tfDnaBoundID = transcriptID + "__" + tf

			# If the tfDnaBoundID is not found in bulkMolecules, do not add a
			# node for this pair - looks like the fitter deletes some of these
			# TF-gene pairs
			try:
				tfDnaBoundIdx = moleculeIDs.index(tfDnaBoundID)
			except ValueError:
				tfDnaBoundIdx = -1

			if tfDnaBoundIdx != -1:
				# Get corresponding ID of gene from the transcript ID
				geneID = gene_ids[rna_ids.index(transcriptID)]

				# Initialize a single regulation node for each TF-gene pair
				regulation_node = Node("Process", "Regulation")

				# Add attributes to the node
				# TODO: Add gene regulation strength constants?
				regID = tf + "_" + geneID + "_REGULATION"
				regName = tf + "-" + geneID + " gene regulation"
				attr = {'node_id': regID,
					'name': regName,
				}
				regulation_node.read_attributes(**attr)

				# Add dynamics data (counts) to the node.
				# Get column index of the TF+gene pair in the counts array
				dynamics = {'bound TFs': list(counts_array[:, tfDnaBoundIdx].astype(np.int))}
				dynamics_units = {'bound TFs': 'N'}
				regulation_node.read_dynamics(dynamics, dynamics_units)

				node_list.append(regulation_node)

				# Add edge from TF to this regulation node
				regulation_edge_from_tf = Edge("Regulation")
				attr = {'src_id': tfID,
					'dst_id': regID,
				}
				regulation_edge_from_tf.read_attributes(**attr)
				edge_list.append(regulation_edge_from_tf)

				# Add edge from this regulation node to the gene
				regulation_edge_to_gene = Edge("Regulation")
				attr = {'src_id': regID,
					'dst_id': geneID,
				}
				regulation_edge_to_gene.read_attributes(**attr)
				edge_list.append(regulation_edge_to_gene)


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

	time_string = format_dynamics_string(time, "time")

	time_row = "%s,%s,%s,%s\n" % (
		"time", "time", "s", time_string
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

		# Format dynamics string depending on data type
		if unit == "N":
			dynamics_string = format_dynamics_string(dynamics, "int")
		else:
			dynamics_string = format_dynamics_string(dynamics, "float")

		# Format single string with dynamic attributes separated by commas
		dynamics_row = "%s,%s,%s,%s\n" % (
			"global", name, unit, dynamics_string
		)

		# Write line to dynamics file
		dynamics_file.write(dynamics_row)


def find_duplicate_nodes(node_list):
	"""
	Identify any nodes that have duplicate IDs and prints them to the console.
	This does not remove duplicate nodes.
	"""
	print("Checking for duplicate node IDs...")

	node_ids = []
	duplicate_ids = []

	# Loop through all nodes in the node_list
	for node in node_list:
		# Get ID of the node
		node_id = node.get_node_id()

		# If node was not seen, add to returned list of unique node IDs
		if node_id not in node_ids:
			node_ids.append(node_id)
		# If node was seen, add to list of duplicate IDs
		elif node_id not in duplicate_ids:
			duplicate_ids.append(node_id)

	# Print duplicate node IDs that were found
	for node_id in duplicate_ids:
		print("Node ID %s is duplicate." % (node_id, ))

	return node_ids


def find_runaway_edges(node_ids, edge_list):
	"""
	Find any edges that connect nonexistent nodes and prints them to the
	console.
	"""
	print("Checking for runaway edges...")

	# Loop through all edges in edge_list
	for edge in edge_list:
		# Get IDs of source node
		src_id = edge.get_src_id()
		dst_id = edge.get_dst_id()
		process = edge.get_process()

		# Print error prompt if the node IDs are not found in node_ids
		if src_id not in node_ids:
			print("src_id %s of a %s edge does not exist." % (src_id, process, ))

		if dst_id not in node_ids:
			print("dst_id %s of a %s edge does not exist." % (dst_id, process, ))


def format_dynamics_string(dynamics, datatype):
	"""
	Formats the string of dynamics data that is printed out to the dynamics
	file. If datatype is "int", print all numbers as full decimal integers.
	If datatype is "float", print all numbers in the general format with the
	precision set by DYNAMICS_PRECISION. If datatype is "time", print all
	numbers in the floating point format with the precision set by
	TIME_PRECISION.
	"""
	if datatype == "int":
		# Format first datapoint
		dynamics_string = "[{0:d}".format(dynamics[0])

		# Format rest of the datapoints
		for val in dynamics[1:]:
			dynamics_string += ", {0:d}".format(val)
		dynamics_string += "]"

	elif datatype == "float":
		# Format first datapoint
		dynamics_string = "[{0:.{1}g}".format(dynamics[0], DYNAMICS_PRECISION)

		# Format rest of the datapoints
		for val in dynamics[1:]:
			dynamics_string += ", {0:.{1}g}".format(val, DYNAMICS_PRECISION)
		dynamics_string += "]"

	elif datatype == "time":
		# Format first datapoint
		dynamics_string = "[{0:.{1}f}".format(dynamics[0], TIME_PRECISION)

		# Format rest of the datapoints
		for val in dynamics[1:]:
			dynamics_string += ", {0:.{1}f}".format(val, TIME_PRECISION)
		dynamics_string += "]"

	else:
		dynamics_string = dynamics

	return dynamics_string


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

	# Check for network sanity (optional)
	if CHECK_SANITY:
		print("Performing sanity check on network...")
		node_ids = find_duplicate_nodes(node_list)
		find_runaway_edges(node_ids, edge_list)
		print("Sanity check completed.")

	print("Total number of nodes: %d" % (len(node_list)))
	print("Total number of edges: %d" % (len(edge_list)))

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
