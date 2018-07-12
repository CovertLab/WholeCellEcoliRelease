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
		self.dynamics = {}
		self.dynamics_units = {}

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
	# Create nodes for dNTPs
	dNTPs_nodelist = []
	for dNTP in simData.moleculeGroups.dNtpIds:
		dNTP_node = Node("State", "Metabolite")
		attr = {'node_id': dNTP,
			'name': dNTP,
			'constants': {'mass': 0},
			}
		dNTP_node.get_attributes(**attr)
		dNTPs_nodelist.append(dNTP_node)

		# Add dNTP nodes to node list
		node_list.append(dNTP_node)

	# Create node for Ppi
	ppi_node = Node("State", "Metabolite")
	attr = {'node_id': 'PPI[c]',
			'name': 'PPI[c]',
			'constants': {'mass': 0},
			}
	ppi_node.get_attributes(**attr)

	# Add Ppi node to node list
	node_list.append(ppi_node)

	# Create nodes for DNA polymerases
	pols_nodelist = []
	for pol in ['CPLX0-2361[c]', 'CPLX0-3761[c]', 'CPLX0-3925[c]', 'CPLX0-7910[c]']:
		pol_node = Node("State", "Complex")
		attr = {'node_id': pol,
			'name': pol,
			'constants': {'mass': 0},
			}
		pol_node.get_attributes(**attr)
		pols_nodelist.append(pol_node)

		# Add DNA polymerase nodes to node list
		node_list.append(pol_node)

	# Loop through all genes
	geneIds = [data[0] for data in simData.process.replication.geneData]
	for geneId in geneIds:
		# Initialize a single gene node
		gene_node = Node("State", "Gene")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': geneId, 'name': geneId}
		gene_node.get_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append gene node to node_list
		node_list.append(gene_node)

		# Initialize a single replication node for each gene
		replication_node = Node("Process", "Replication")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': "%s_REPLICATION" % geneId, 'name': "%s_REPLICATION" % geneId}
		replication_node.get_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append replication node to node_list
		node_list.append(replication_node)

		# Add edges between gene and replication nodes
		gene_to_replication_edge = Edge("Replication")
		attr = {'src_id': gene_node, 'dst_id': replication_node}
		gene_to_replication_edge.get_attributes(**attr)
		edge_list.append(gene_to_replication_edge)

		replication_to_gene_edge = Edge("Replication")
		attr = {'src_id': replication_node, 'dst_id': gene_node}
		gene_to_replication_edge.get_attributes(**attr)
		edge_list.append(replication_to_gene_edge)

		# Add edges from dNTPs to replication nodes
		for dNTP_node in dNTPs_nodelist:
			dNTP_to_replication_edge = Edge("Replication")
			attr = {'src_id': dNTP_node, 'dst_id': replication_node}
			dNTP_to_replication_edge.get_attributes(**attr)
			edge_list.append(dNTP_to_replication_edge)

		# Add edge from replication to Ppi
		replication_to_ppi_edge = Edge("Replication")
		attr = {'src_id': replication_node, 'dst_id': ppi_node}
		replication_to_ppi_edge.get_attributes(**attr)
		edge_list.append(replication_to_ppi_edge)

		# Add edges from DNA polymerases to replication
		for pol_node in pols_nodelist:
			pol_to_replication_edge = Edge("Replication")
			attr = {'src_id': pol_node, 'dst_id': replication_node}
			pol_to_replication_edge.get_attributes(**attr)
			edge_list.append(pol_to_replication_edge)


def add_transcription_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add transcription nodes with dynamics data to the node list, and add edges
	connected to the transcription nodes to the edge list. - Heejo
	"""
	# Create nodes for NTPs
	ntps_nodelist = []
	for ntp in simData.moleculeGroups.ntpIds:
		ntp_node = Node("State", "Metabolite")
		attr = {'node_id': ntp,
			'name': ntp,
			'constants': {'mass': 0},
			}
		ntp_node.get_attributes(**attr)
		ntps_nodelist.append(ntp_node)

		# Add NTP node to node list
		node_list.append(ntp_node)

	# Create node for Ppi
	ppi_node = Node("State", "Metabolite")
	attr = {'node_id': 'PPI[c]',
			'name': 'PPI[c]',
			'constants': {'mass': 0},
			}
	ppi_node.get_attributes(**attr)

	# Add Ppi node to node list
	node_list.append(ppi_node)

	# Create node for RNA polymerase
	pol_node = Node("State", "Complex")
	attr = {'node_id': 'APORNAP-CPLX[c]',
		'name': 'APORNAP-CPLX[c]',
		'constants': {'mass': 0},
		}
	pol_node.get_attributes(**attr)

	# Add RNA polymerase nodes to node list
	node_list.append(pol_node)

	# Loop through all genes
	for geneId, rnaId, _ in simData.process.replication.geneData:
		# Initialize a single gene node
		gene_node = Node("State", "Gene")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': geneId, 'name': geneId}
		gene_node.get_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append gene node to node_list
		node_list.append(gene_node)

		# Initialize a single transcript node
		rna_node = Node("State", "RNA")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': "%s[c]" % rnaId, 'name': "%s[c]" % rnaId}
		rna_node.get_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append transcript node to node_list
		node_list.append(rna_node)

		# Initialize a single transcription node for each gene
		transcription_node = Node("Process", "Transcription")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': "%s_TRANSCRIPTION" % geneId, 'name': "%s_TRANSCRIPTION" % geneId}
		transcription_node.get_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append transcription node to node_list
		node_list.append(transcription_node)

		# Add edge from gene to transcription node
		gene_to_transcription_edge = Edge("Transcription")
		attr = {'src_id': gene_node, 'dst_id': transcription_node}
		gene_to_transcription_edge.get_attributes(**attr)
		edge_list.append(gene_to_transcription_edge)

		# Add edge from transcription to transcript node
		transcription_to_rna_edge = Edge("Transcription")
		attr = {'src_id': transcription_node, 'dst_id': rna_node}
		transcription_to_rna_edge.get_attributes(**attr)
		edge_list.append(transcription_to_rna_edge)

		# Add edges from NTPs to transcription nodes
		for ntp_node in ntps_nodelist:
			ntp_to_transcription_edge = Edge("Transcription")
			attr = {'src_id': ntp_node, 'dst_id': transcription_node}
			ntp_to_transcription_edge.get_attributes(**attr)
			edge_list.append(ntp_to_transcription_edge)

		# Add edge from transcription to Ppi
		transcription_to_ppi_edge = Edge("Transcription")
		attr = {'src_id': transcription_node, 'dst_id': ppi_node}
		transcription_to_ppi_edge.get_attributes(**attr)
		edge_list.append(transcription_to_ppi_edge)

		# Add edges from RNA polymerases to transcription
		pol_to_transcription_edge = Edge("Transcription")
		attr = {'src_id': pol_node, 'dst_id': transcription_node}
		pol_to_transcription_edge.get_attributes(**attr)
		edge_list.append(pol_to_transcription_edge)


def add_translation_nodes_and_edges(simData, simOutDirs, node_list, edge_list):
	"""
	Add translation nodes with dynamics data to the node list, and add edges
	connected to the translation nodes to the edge list. - Heejo
	"""
	# Create nodes for amino acids
	aas_nodelist = []
	for aa in simData.moleculeGroups.aaIDs:
		aa_node = Node("State", "Metabolite")
		attr = {'node_id': aa,
			'name': aa,
			'constants': {'mass': 0},
			}
		aa_node.get_attributes(**attr)
		aas_nodelist.append(aa_node)

		# Add amino acid node to node list
		node_list.append(aa_node)

	# Create nodes for GTP, GDP, water, Ppi
	metabolites_nodelist = []
	for metabolite in ["GTP[c]", "GDP[c]", "WATER[c]", "PPI[c]"]:
		metabolite_node = Node("State", "Metabolite")
		attr = {'node_id': metabolite,
			'name': metabolite,
			'constants': {'mass': 0},
			}
		metabolite_node.get_attributes(**attr)
		metabolites_nodelist.append(metabolite_node)

		# Add metabolite node to node list
		node_list.append(metabolite_node)
	gtp_node, gdp_node, water_node, ppi_node = metabolites_nodelist

	# Create nodes for ribosome subunits
	ribosomes_nodelist = []
	subunits = [simData.moleculeGroups.s30_fullComplex, simData.moleculeGroups.s50_fullComplex]
	for subunit in subunits:
		subunit_node = Node("State", "Complex")
		attr = {'node_id': subunit,
			'name': subunit,
			'constants': {'mass': 0},
			}
		subunit_node.get_attributes(**attr)
		ribosomes_nodelist.append(subunit_node)

		# Add ribosome subunits to node list
		node_list.append(subunit_node)

	# Loop through all translatable genes
	for data in simData.process.translation.monomerData:
		monomerId = data[0]
		rnaId = data[1]
		geneId = rnaId.split("_RNA[c]")[0]

		# Initialize a single transcript node
		rna_node = Node("State", "RNA")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': rnaId, 'name': rnaId} # rnaId already contains [c] tag
		rna_node.get_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append transcript node to node_list
		node_list.append(rna_node)

		# Initialize a single protein node
		protein_node = Node("State", "Protein")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': monomerId, 'name': monomerId}
		protein_node.get_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append protein node to node_list
		node_list.append(protein_node)

		# Initialize a single translation node for each transcript
		translation_node = Node("Process", "Translation")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': "%s_TRANSLATION" % geneId, 'name': "%s_TRANSLATION" % geneId}
		translation_node.get_attributes(**attr)

		# Add dynamics data to the node.
		# TODO

		# Append translation node to node_list
		node_list.append(translation_node)

		# Add edge from transcript to translation node
		rna_to_translation_edge = Edge("Translation")
		attr = {'src_id': rna_node, 'dst_id': translation_node}
		rna_to_translation_edge.get_attributes(**attr)
		edge_list.append(rna_to_translation_edge)

		# Add edge from translation to monomer node
		translation_to_protein_edge = Edge("Translation")
		attr = {'src_id': translation_node, 'dst_id': protein_node}
		translation_to_protein_edge.get_attributes(**attr)
		edge_list.append(translation_to_protein_edge)

		# Add edges from amino acids to translation node
		for aa_node in aas_nodelist:
			aa_to_translation_edge = Edge("Translation")
			attr = {'src_id': aa_node, 'dst_id': translation_node}
			aa_to_translation_edge.get_attributes(**attr)
			edge_list.append(aa_to_translation_edge)

		# Add edges from other reactants to translation node
		for reactant_node in [gtp_node, water_node]:
			reactant_to_translation_edge = Edge("Translation")
			attr = {'src_id': reactant_node, 'dst_id': translation_node}
			reactant_to_translation_edge.get_attributes(**attr)
			edge_list.append(reactant_to_translation_edge)

		# Add edges from translation to other product nodes
		for product_node in [gdp_node, ppi_node, water_node]:
			translation_to_product_edge = Edge("Translation")
			attr = {'src_id': translation_node, 'dst_id': product_node}
			translation_to_product_edge.get_attributes(**attr)
			edge_list.append(translation_to_product_edge)

		# Add edges from ribosome subunits to translation node
		for subunit_node in ribosomes_nodelist:
			subunit_to_translation_edge = Edge("Translation")
			attr = {'src_id': subunit_node, 'dst_id': translation_node}
			subunit_to_translation_edge.get_attributes(**attr)
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

	# Initialize list of protein and complex IDs
	protein_ids = []
	complex_ids = []

	# Loop through all complexation reactions
	for idx, reaction in enumerate(reactionIDs):
		# Initialize a single complexation node for each complexation reaction
		complexation_node = Node("Process", "Complexation")

		# Add attributes to the node
		attr = {'node_id': reaction, 'name': reaction}
		complexation_node.get_attributes(**attr)

		# # TODO (ERAN) Add dynamics data (# rxns/ts) to the node.
		# dynamics = {'flux': list(flux_array[:, idx])}
		# dynamics_units = {'flux': 'mmol/gCDW/h'}
		# complexation_node.get_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(complexation_node)

		# Get reaction stoichiometry from complexStoich
		stoich_dict = complexStoich[reaction]

		# Loop through all proteins participating in the reaction
		for protein, stoich in stoich_dict.items():
			# Add complexes that were not encountered
			if ('CPLX' in protein) and (protein not in complex_ids):
				complex_ids.append(protein)
			# Add proteins that were not encountered
			elif protein not in protein_ids:
				protein_ids.append(protein)

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
			complex_edge.get_attributes(**attr)

			# Append edge to edge_list
			edge_list.append(complex_edge)

	# Loop through all proteins
	for protein in protein_ids:
		# Initialize a single protein node for each protein
		protein_node = Node("State", "Protein")

		# Add attributes to the node
		# TODO: Get molecular mass using getMass().
		# TODO: Get correct protein name and synonyms from EcoCyc
		attr = {'node_id': protein,
			'name': protein,
			'constants': {'mass': 0}
			}
		protein_node.get_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the protein in the counts array
		try:
			protein_idx = moleculeIDs.index(protein)
		except ValueError:  # protein ID not found in moleculeIDs
			protein_idx = -1

		if protein_idx != -1:
			dynamics = {'counts': list(counts_array[:, protein_idx])}
			dynamics_units = {'counts': 'N'}
			protein_node.get_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(protein_node)


	# Loop through all complexes
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
		complex_node.get_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the complex in the counts array
		try:
			complex_idx = moleculeIDs.index(complex)
		except ValueError:  # complex ID not found in moleculeIDs
			complex_idx = -1

		if complex_idx != -1:
			dynamics = {'counts': list(counts_array[:, complex_idx])}
			dynamics_units = {'counts': 'N'}
			complex_node.get_dynamics(dynamics, dynamics_units)

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
	counts_array = np.empty((0, len(moleculeIDs)))

	for simOutDir in simOutDirs:
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionFluxes = fbaResults.readColumn('reactionFluxes')
		flux_array = np.concatenate((flux_array, reactionFluxes))

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

	# Get all reaction stoichiometry from simData
	reactionStoich = simData.process.metabolism.reactionStoich

	# Initialize list of metabolite IDs
	metabolite_ids = []

	# Loop through all reactions
	for idx, reaction in enumerate(reactionIDs):
		# Initialize a single metabolism node for each reaction
		metabolism_node = Node("Process", "Metabolism")

		# Add attributes to the node
		# TODO: Get correct reaction name and synonyms from EcoCyc
		attr = {'node_id': reaction, 'name': reaction}
		metabolism_node.get_attributes(**attr)

		# Add dynamics data (flux) to the node. The flux array shares the same
		# column index with the reactionIDs list
		dynamics = {'flux': list(flux_array[:, idx])}
		dynamics_units = {'flux': 'mmol/gCDW/h'}
		metabolism_node.get_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(metabolism_node)

		# Get reaction stoichiometry from reactionStoich
		stoich_dict = reactionStoich[reaction]

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
			metabolism_edge.get_attributes(**attr)

			# Append edge to edge_list
			edge_list.append(metabolism_edge)

	# Loop through all metabolites
	for metabolite in metabolite_ids:
		# Initialize a single metabolite node for each metabolite
		metabolite_node = Node("State", "Metabolite")

		# Add attributes to the node
		# TODO: Get molecular mass using getMass(). Some of the metabolites do not have mass data?
		# TODO: Get correct metabolite name and synonyms from EcoCyc
		attr = {'node_id': metabolite,
			'name': metabolite,
			'constants': {'mass': 0}
			}
		metabolite_node.get_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the metabolite in the counts array
		try:
			metabolite_idx = moleculeIDs.index(metabolite)
		except ValueError:  # metabolite ID not found in moleculeIDs
			# Some of the metabolites are not being tracked.
			# TODO: why?
			metabolite_idx = -1

		if metabolite_idx != -1:
			dynamics = {'counts': list(counts_array[:, metabolite_idx])}
			dynamics_units = {'counts': 'N'}
			metabolite_node.get_dynamics(dynamics, dynamics_units)

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

	# TODO: get dynamics data for equlibrium nodes. Will need new listener.

	# Loop through each equilibrium reaction
	for reactionIdx, rxnId in enumerate(rxnIds):
		# Initialize a single equilibrium node for each equilibrium reaction
		equilibrium_node = Node("Process", "Equilibrium")

		# Add attributes to the node
		# TODO: Think of better "name" for this node?
		attr = {'node_id': rxnId,
			'name': rxnId,
			'constants': {'rateFwd': ratesFwd[reactionIdx], 'rateRev': ratesRev[reactionIdx]}
		}
		equilibrium_node.get_attributes(**attr)

		# Append new node to node_list
		node_list.append(equilibrium_node)

		# Extract column corresponding to reaction in the stoichiometric matrix
		stoichMatrixColumn = stoichMatrix[:, reactionIdx]

		# Loop through each element in column
		for moleculeIdx, stoich in enumerate(stoichMatrixColumn):
			# If the stoichiometric coefficient is negative, add reactant edge
			# to the equilibrium node
			if stoich < 0:
				equilibrium_edge = Edge("Equilibrium")
				attr = {'src_id': moleculeIds[moleculeIdx],
					'dst_id': rxnId,
					'stoichiometry': stoich
				}

				equilibrium_edge.get_attributes(**attr)
				edge_list.append(equilibrium_edge)

			# If the coefficient is positive, add product edge
			elif stoich > 0:
				equilibrium_edge = Edge("Equilibrium")
				attr = {'src_id': rxnId,
					'dst_id': moleculeIds[moleculeIdx],
					'stoichiometry': stoich
				}

				equilibrium_edge.get_attributes(**attr)
				edge_list.append(equilibrium_edge)


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
	# check if all IDs in dynamics correspond to nodes.
	# check if all nodes have at least one edge.
	# check that all edges connect nodes that are in the nodelist.
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
