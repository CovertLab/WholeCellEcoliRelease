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
import csv

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

NODE_LIST_HEADER = "ID\tclass\tcategory\tname\tsynonyms\tconstants\n"
EDGE_LIST_HEADER = "src_node_id\tdst_node_id\tstoichiometry\tprocess\n"
DYNAMICS_HEADER = "node\ttype\tunits\tdynamics\n"
PATHWAY_LIST_HEADER = "pathway\tnodes"

PATHWAYS_FILENAME = "models/ecoli/analysis/causal_network/metabolic_pathways.tsv"

NAMES_PATHWAY = "models/ecoli/analysis/causal_network/names/"

CHECK_SANITY = False
GET_PATHWAY_INDEX = False
N_GENS = 9 # TODO (Eran) this is structures as a multigen analysis, what if we want to analyze single gen?
DYNAMICS_PRECISION = 6
PROBABILITY_PRECISION = 4
TIME_PRECISION = 2

# Proteins that are reactants and products of a metabolic reaction
PROTEINS_IN_METABOLISM = ["EG50003-MONOMER[c]", "PHOB-MONOMER[c]", "PTSI-MONOMER[c]", "PTSH-MONOMER[c]"]

# Equilibrium complexes that are formed from deleted equilibrium reactions, but
# are reactants in a complexation reaction
EQUILIBRIUM_COMPLEXES_IN_COMPLEXATION = ["CPLX0-7620[c]", "CPLX0-7701[c]", "CPLX0-7677[c]", "MONOMER0-1781[c]", "CPLX0-7702[c]"]

# Metabolites that are used as ligands in equilibrium, but do not participate
# in any metabolic reactions
METABOLITES_ONLY_IN_EQUILIBRIUM = ["4FE-4S[c]", "NITRATE[p]"]

# Molecules in 2CS reactions that are not proteins
NONPROTEIN_MOLECULES_IN_2CS = ["ATP[c]", "ADP[c]", "WATER[c]", "PI[c]", "PROTON[c]", "PHOSPHO-PHOB[c]"]

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
		self.names_dict = {}

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
				dynamics_string = format_dynamics_string(dynamics, "int")
			elif unit == "prob":
				dynamics_string = format_dynamics_string(dynamics, "prob")
			else:
				dynamics_string = format_dynamics_string(dynamics, "float")

			# Format single string with dynamic attributes separated by commas
			dynamics_row = "%s\t%s\t%s\t%s\n" % (
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
		edge_row = "%s\t%s\t%s\t%s\n" % (
			self.src_id, self.dst_id, self.stoichiometry, self.process,
			)

		# Write line to edgelist file
		edgelist_file.write(edge_row)


def add_global_nodes(simData, simOutDirs, node_list):
	"""
	Add global state nodes to the node list.
	"""
	# Get bulkMolecule IDs from first simOut directory
	simOutDir = simOutDirs[0]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get index of full chromosome
	full_chrom_idx = moleculeIDs.index(simData.moleculeGroups.fullChromosome[0])

	mass_array = np.empty(0)
	volume_array = np.empty(0)
	fc_counts_array = np.empty(0, dtype=np.int)

	# Loop through all generations
	for simOutDir in simOutDirs:
		# Extract dynamics data from each generation
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")

		cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
		cell_volume = ((1.0/simData.constants.cellDensity)*(units.fg*cell_mass)).asNumber(units.L)

		# Append to existing array
		fc_counts_array = np.concatenate((fc_counts_array, counts[:, full_chrom_idx].astype(np.int)))
		mass_array = np.concatenate((mass_array, cell_mass))
		volume_array = np.concatenate((volume_array, cell_volume))

	# Add total cell mass node to node list
	mass_node = Node("State", "Global")
	attr = {'node_id': "global", 'name': "Total cell mass"}
	mass_node.read_attributes(**attr)

	dynamics = {'mass': list(mass_array)}
	dynamics_units = {'mass': 'fg'}
	mass_node.read_dynamics(dynamics, dynamics_units)

	# Add total cell volume node to node list
	volume_node = Node("State", "Global")
	attr = {'node_id': "global", 'name': "Total cell volume"}
	volume_node.read_attributes(**attr)

	dynamics = {'volume': list(volume_array)}
	dynamics_units = {'volume': 'L'}
	volume_node.read_dynamics(dynamics, dynamics_units)

	# Add chromosome count node to node list
	chromosome_node = Node("State", "Global")
	attr = {'node_id': "global", 'name': "Full chromosome counts"}
	chromosome_node.read_attributes(**attr)

	dynamics = {'count': list(fc_counts_array)}
	dynamics_units = {'count': 'N'}
	chromosome_node.read_dynamics(dynamics, dynamics_units)

	node_list.extend([mass_node, volume_node, chromosome_node])


def add_replication_and_genes(simData, simOutDirs, node_list, edge_list, names_dict):
	"""
	Add replication process nodes and gene state nodes with dynamics data to
	the node list, and the edges connected to the replication nodes to the edge
	list. - Heejo
	"""
	dntp_ids = simData.moleculeGroups.dNtpIds
	ppi_id = "PPI[c]"
	dnap_ids = ['CPLX0-2361[c]', 'CPLX0-3761[c]', 'CPLX0-3925[c]', 'CPLX0-7910[c]']

	# Assemble dynamics data from all generations
	rnaSynthProb = []
	for simOutDir in simOutDirs:
		rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		rnaSynthProb_thisGen = rnaSynthProbReader.readColumn("rnaSynthProb")
		rnaSynthProbReader.close()
		rnaSynthProb.append(rnaSynthProb_thisGen)
	rnaSynthProb = np.concatenate(rnaSynthProb)

	# Loop through all genes (in the order listed in transcription)
	for i, geneId in enumerate(simData.process.transcription.rnaData["geneId"]):
		# Initialize a single gene node
		gene_node = Node("State", "Gene")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': geneId, 'name': geneId}
		gene_node.read_attributes(**attr)

		# Add dynamics data to the node. The rna synthesis probability shares
		# the same index as the geneId.
		dynamics = {"transcript synthesis probability": list(rnaSynthProb[:, i])}
		dynamics_units = {"transcript synthesis probability": "p"}
		gene_node.read_dynamics(dynamics, dynamics_units)

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


def add_transcription_and_transcripts(simData, simOutDirs, node_list, edge_list, names_dict):
	"""
	Add transcription process nodes and transcript state nodes with dynamics
	data to the node list, and edges connected to the transcription nodes to
	the edge list. - Heejo
	"""
	ntp_ids = simData.moleculeGroups.ntpIds
	ppi_id = "PPI[c]"
	rnap_id = "APORNAP-CPLX[c]"

	# Get bulkMolecule IDs from first simOut directory
	simOutDir = simOutDirs[0]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get dynamics data from all simOutDirs
	counts_array = np.empty((0, len(moleculeIDs)), dtype=np.int)
	nRnaInits = []

	for simOutDir in simOutDirs:
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

		rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
		rnaInitEvent = rnapDataReader.readColumn("rnaInitEvent")
		rnapDataReader.close()
		nRnaInits.append(rnaInitEvent)
	nRnaInits = np.concatenate(nRnaInits)

	# Loop through all genes (in the order listed in transcription)
	for i, rnaId in enumerate(simData.process.transcription.rnaData["id"]):
		geneId = simData.process.transcription.rnaData["geneId"][i]

		# Initialize a single transcript node
		rna_node = Node("State", "RNA")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': rnaId, 'name': rnaId}
		rna_node.read_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the RNA in the counts array
		try:
			rna_idx = moleculeIDs.index(rnaId)
		except ValueError:  # RNA ID not found in moleculeIDs
			rna_idx = -1

		if rna_idx != -1:
			rna_counts = counts_array[:, rna_idx].astype(np.int)

			dynamics = {'counts': list(rna_counts)}
			dynamics_units = {'counts': 'N'}
			rna_node.read_dynamics(dynamics, dynamics_units)

		# Append transcript node to node_list
		node_list.append(rna_node)

		# Initialize a single transcription node for each gene
		transcription_node = Node("Process", "Transcription")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		transcription_node_id = "%s_TRANSCRIPTION" % geneId
		attr = {'node_id': transcription_node_id, 'name': transcription_node_id}
		transcription_node.read_attributes(**attr)

		# Add dynamics data to the node. The number of transcription initiation
		# events (per gene per second) shares the same index as the rnaId.
		dynamics = {"transcription initiations": list(nRnaInits[:, i])}
		dynamics_units = {"transcription initiations": "N"}
		transcription_node.read_dynamics(dynamics, dynamics_units)

		# Append transcription node to node_list
		node_list.append(transcription_node)

		# Add edge from gene to transcription node
		gene_to_transcription_edge = Edge("Transcription")
		attr = {'src_id': geneId, 'dst_id': transcription_node_id}
		gene_to_transcription_edge.read_attributes(**attr)
		edge_list.append(gene_to_transcription_edge)

		# Add edge from transcription to transcript node
		transcription_to_rna_edge = Edge("Transcription")
		attr = {'src_id': transcription_node_id, 'dst_id': rnaId}
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


def add_translation_and_monomers(simData, simOutDirs, node_list, edge_list, names_dict):
	"""
	Add translation process nodes and protein (monomer) state nodes with
	dynamics data to the node list, and edges connected to the translation
	nodes to the edge list. - Heejo
	"""
	# Create nodes for amino acids
	aa_ids = simData.moleculeGroups.aaIDs
	gtp_id = "GTP[c]"
	gdp_id = "GDP[c]"
	water_id = "WATER[c]"
	ppi_id = "PPI[c]"

	ribosome_subunit_ids = [simData.moleculeGroups.s30_fullComplex[0], simData.moleculeGroups.s50_fullComplex[0]]
	n_avogadro = simData.constants.nAvogadro

	# Get bulkMolecule IDs from first simOut directory
	simOutDir = simOutDirs[0]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	nMonomers = len(simData.process.translation.monomerData)

	# Get dynamics data from all simOutDirs
	probTranslation_array = np.empty((0, nMonomers), dtype=np.int)
	counts_array = np.empty((0, len(moleculeIDs)), dtype=np.int)
	volume_array = np.empty(0)

	for simOutDir in simOutDirs:
		translationResults = TableReader(os.path.join(simOutDir, "RibosomeData"))
		probTranslationPerTranscript = translationResults.readColumn('probTranslationPerTranscript')
		probTranslation_array = np.concatenate((probTranslation_array, probTranslationPerTranscript))

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

		# Extract dynamics data from each generation
		cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
		cell_volume = ((1.0/simData.constants.cellDensity)*(units.fg*cell_mass)).asNumber(units.L)
		volume_array = np.concatenate((volume_array, cell_volume))

	# Loop through all translatable genes
	for idx, data in enumerate(simData.process.translation.monomerData):
		monomerId = data[0]
		rnaId = data[1]
		geneId = rnaId.split("_RNA[c]")[0]

		# Initialize a single protein node
		protein_node = Node("State", "Protein")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		attr = {'node_id': monomerId, 'name': monomerId}
		protein_node.read_attributes(**attr)

		# Add dynamics data (counts) to the node.
		# Get column index of the monomer in the counts array
		try:
			monomer_idx = moleculeIDs.index(monomerId)
		except ValueError:  # complex ID not found in moleculeIDs
			monomer_idx = -1

		if monomer_idx != -1:
			monomer_counts = counts_array[:, monomer_idx].astype(np.int)
			monomer_conc = (((1/n_avogadro)*monomer_counts)/(units.L*volume_array)).asNumber(units.mmol/units.L)

			dynamics = {'counts': list(monomer_counts),
				'concentration': list(monomer_conc)}
			dynamics_units = {'counts': 'N', 'concentration': 'mmol/L'}
			protein_node.read_dynamics(dynamics, dynamics_units)

		# Append protein node to node_list
		node_list.append(protein_node)

		# Initialize a single translation node for each transcript
		translation_node = Node("Process", "Translation")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		translation_node_id = "%s_TRANSLATION" % geneId
		attr = {'node_id': translation_node_id, 'name': translation_node_id}
		translation_node.read_attributes(**attr)

		# Add dynamics data (probability of translation initiation per transcript) to the node.
		dynamics = {'probability of initiating translation': list(probTranslation_array[:, idx])}
		dynamics_units = {'probability of initiating translation': 'p'}
		translation_node.read_dynamics(dynamics, dynamics_units)

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


def add_complexation_and_complexes(simData, simOutDirs, node_list, edge_list, names_dict):
	"""
	Add complexation process nodes and complex state nodes with dynamics data
	to the node list, and edges connected to the complexation nodes to the edge
	list. - Eran
	"""
	simOutDir = simOutDirs[0]

	# TODO (Eran) raw_data is here used to get complexation stoichiometry. This can be saved to sim_data and then retrieved here.
	from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
	raw_data = KnowledgeBaseEcoli()

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	n_avogadro = simData.constants.nAvogadro

	# List of all complex IDs and reaction IDs
	complex_ids = simData.process.complexation.ids_complexes + EQUILIBRIUM_COMPLEXES_IN_COMPLEXATION
	reactionIDs = simData.process.complexation.ids_reactions

	# Get dynamics data from all simOutDirs (# rxns/ts for complexation, counts)
	reactions_array = np.empty((0, len(reactionIDs)))
	counts_array = np.empty((0, len(moleculeIDs)))
	volume_array = np.empty(0)

	for simOutDir in simOutDirs:
		complexationResults = TableReader(os.path.join(simOutDir, "ComplexationListener"))
		reactionRates = complexationResults.readColumn('reactionRates')
		reactions_array = np.concatenate((reactions_array, reactionRates))

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

		# Extract dynamics data from each generation
		cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
		cell_volume = ((1.0 / simData.constants.cellDensity) * (units.fg * cell_mass)).asNumber(units.L)
		volume_array = np.concatenate((volume_array, cell_volume))

	# Get complexation stoichiometry from simData
	complexStoich = {}
	# TODO (Eran) save complexationReactions in sim_data, so that raw_data won't be needed
	for reaction in raw_data.complexationReactions:
		stoich = {}
		for molecule in reaction['stoichiometry']:
			molecule_name = '%s[%s]' % (molecule['molecule'], molecule['location'])
			stoich[molecule_name] = molecule['coeff']
		complexStoich[reaction['id']] = stoich

	# Loop through all complexation reactions
	for idx, reaction in enumerate(reactionIDs):
		# Initialize a single complexation node for each complexation reaction
		complexation_node = Node("Process", "Complexation")

		# Add attributes to the node
		attr = {'node_id': reaction, 'name': reaction}
		complexation_node.read_attributes(**attr)

		# Add dynamics data (# rxns/sec) to the node.
		dynamics = {'reaction rate': list(reactions_array[:, idx])}
		dynamics_units = {'reaction rate': 'rxns/s'}
		complexation_node.read_dynamics(dynamics, dynamics_units)

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
			complex_counts = counts_array[:, complex_idx].astype(np.int)
			complex_conc = (((1 / n_avogadro) * complex_counts) / (units.L * volume_array)).asNumber(units.mmol / units.L)

			dynamics = {'counts': list(complex_counts),
				'concentration':  list(complex_conc)}
			dynamics_units = {'counts': 'N', 'concentration': 'mmol/L'}
			complex_node.read_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(complex_node)


def add_metabolism_and_metabolites(simData, simOutDirs, node_list, edge_list, names_dict):
	"""
	Add metabolism process nodes and metabolite state nodes with dynamics data
	to the node list, add edges connected to the metabolism nodes to the edge
	list. - Gwanggyu
	Note: forward and reverse reactions are represented as separate nodes.
	"""
	# Get reaction list from first simOut directory
	simOutDir = simOutDirs[0]
	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = fbaResults.readAttribute("reactionIDs")

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	n_avogadro = simData.constants.nAvogadro

	# Get dynamics data from all simOutDirs (flux for reactions, counts for
	# metabolites)
	flux_array = np.empty((0, len(reactionIDs)))
	counts_array = np.empty((0, len(moleculeIDs)), dtype=np.int)
	volume_array = np.empty(0)

	for simOutDir in simOutDirs:
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		reactionFluxes = fbaResults.readColumn('reactionFluxes')
		flux_array = np.concatenate((flux_array, reactionFluxes))

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

		# Extract dynamics data from each generation
		cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
		cell_volume = ((1.0 / simData.constants.cellDensity) * (units.fg * cell_mass)).asNumber(units.L)
		volume_array = np.concatenate((volume_array, cell_volume))

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
			metabolite_counts = counts_array[:, metabolite_idx].astype(np.int)
			metabolite_conc = (((1 / n_avogadro) * metabolite_counts) / (units.L * volume_array)).asNumber(units.mmol / units.L)

			dynamics = {'counts': list(metabolite_counts),
				'concentration':  list(metabolite_conc)}
			dynamics_units = {'counts': 'N', 'concentration': 'mmol/L'}
			metabolite_node.read_dynamics(dynamics, dynamics_units)

		# Append node to node_list
		node_list.append(metabolite_node)


def add_equilibrium(simData, simOutDirs, node_list, edge_list, names_dict):
	"""
	Add equilibrium nodes with dynamics data to the node list, and add edges
	connected to the equilibrium nodes to the edge list. - Gwanggyu
	"""
	# Get equilibrium-specific data from simData
	equilibriumMoleculeIds = simData.process.equilibrium.moleculeNames
	equilibriumRxnIds = simData.process.equilibrium.rxnIds
	equilibriumStoichMatrix = simData.process.equilibrium.stoichMatrix()
	equilibriumRatesFwd = np.array(simData.process.equilibrium.ratesFwd, dtype=np.float32)
	equilibriumRatesRev = np.array(simData.process.equilibrium.ratesRev, dtype=np.float32)

	# Get transcription factor-specific data from simData
	recruitmentColNames = simData.process.transcription_regulation.recruitmentColNames
	tf_ids = sorted(set([x.split("__")[-1] for x in recruitmentColNames if
		x.split("__")[-1] != "alpha"]))
	tfToTfType = simData.process.transcription_regulation.tfToTfType

	# Get bulkMolecule IDs from first simOut directory
	simOutDir = simOutDirs[0]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get dynamics data from all simOutDirs
	counts_array = np.empty((0, len(moleculeIDs)), dtype=np.int)
	pPromoterBoundArray = np.empty((0, len(tf_ids)))
	reactions_array = np.empty((0, len(equilibriumRxnIds)))

	for simOutDir in simOutDirs:
		equilibriumResults = TableReader(os.path.join(simOutDir, "EquilibriumListener"))
		reactionRates = equilibriumResults.readColumn('reactionRates')
		reactions_array = np.concatenate((reactions_array, reactionRates))

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		counts = bulkMolecules.readColumn("counts")
		counts_array = np.concatenate((counts_array, counts))

		rnaSynthProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		pPromoterBound = rnaSynthProb.readColumn("pPromoterBound")
		pPromoterBoundArray = np.concatenate((pPromoterBoundArray, pPromoterBound))

	# Get IDs of complexes that were already added
	complexation_complex_ids = simData.process.complexation.ids_complexes

	# Get list of complex IDs in equilibrium
	equilibrium_complex_ids = simData.process.equilibrium.ids_complexes

	# Loop through each equilibrium reaction
	for reactionIdx, rxnId in enumerate(equilibriumRxnIds):

		# Initialize a single equilibrium node for each equilibrium reaction
		equilibrium_node = Node("Process", "Equilibrium")

		# Add attributes to the node
		rxnName = rxnId[:-4] + " equilibrium reaction"
		attr = {'node_id': rxnId,
			'name': rxnName,
			'constants': {'rateFwd': equilibriumRatesFwd[reactionIdx], 'rateRev': equilibriumRatesRev[reactionIdx]}
		}
		equilibrium_node.read_attributes(**attr)

		# Add dynamics data (# rxns/sec) to the node.
		dynamics = {'reaction rate': list(reactions_array[:, reactionIdx])}
		dynamics_units = {'reaction rate': 'rxns/s'}
		equilibrium_node.read_dynamics(dynamics, dynamics_units)

		# Append new node to node_list
		node_list.append(equilibrium_node)

		# Extract column corresponding to reaction in the stoichiometric matrix
		equilibriumStoichMatrixColumn = equilibriumStoichMatrix[:, reactionIdx]

		# Loop through each element in column
		for moleculeIdx, stoich in enumerate(equilibriumStoichMatrixColumn):
			moleculeId = equilibriumMoleculeIds[moleculeIdx]

			if moleculeId[:-3] in tf_ids and tfToTfType[moleculeId[:-3]] != '0CS':
				tf_idx = tf_ids.index(moleculeId[:-3])

				dynamics['fraction active TF'] = list(pPromoterBoundArray[:, tf_idx])
				dynamics_units['fraction active TF'] = 'prob'

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

		equilibrium_node.read_dynamics(dynamics, dynamics_units)

	# Get 2CS-specific data from simData
	tcsMoleculeIds = simData.process.two_component_system.moleculeNames
	tcsRxnIds = simData.process.two_component_system.rxnIds
	tcsStoichMatrix = simData.process.two_component_system.stoichMatrix()
	tcsRatesFwd = np.array(simData.process.two_component_system.ratesFwd, dtype=np.float32)
	tcsRatesRev = np.array(simData.process.two_component_system.ratesRev, dtype=np.float32)

	# Initialize list of complex IDs in 2CS (should need instance variable)
	tcs_complex_ids = []

	# Get lists of monomers that were already added
	monomer_ids = []
	for monomerData in simData.process.translation.monomerData:
		monomer_ids.append(monomerData[0])

	# Loop through each 2CS reaction
	for reactionIdx, rxnId in enumerate(tcsRxnIds):

		# Initialize a single equilibrium node for each equilibrium reaction
		equilibrium_node = Node("Process", "Equilibrium")

		# Add attributes to the node
		rxnName = rxnId[:-4] + " two-component system reaction"
		attr = {'node_id': rxnId,
			'name': rxnName,
			'constants': {'rateFwd': tcsRatesFwd[reactionIdx], 'rateRev': tcsRatesRev[reactionIdx]}
		}
		equilibrium_node.read_attributes(**attr)

		# Append new node to node_list
		node_list.append(equilibrium_node)

		# Extract column corresponding to reaction in the stoichiometric matrix
		tcsStoichMatrixColumn = tcsStoichMatrix[:, reactionIdx]

		# Loop through each element in column
		for moleculeIdx, stoich in enumerate(tcsStoichMatrixColumn):
			moleculeId = tcsMoleculeIds[moleculeIdx]

			if moleculeId not in monomer_ids + NONPROTEIN_MOLECULES_IN_2CS:
				tcs_complex_ids.append(moleculeId)

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

	# Add new complexes that were encountered here
	for complex_id in list(set(equilibrium_complex_ids + tcs_complex_ids)):
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

	# Loop through metabolites that only appear in equilibrium
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


def add_regulation(simData, simOutDirs, node_list, edge_list, names_dict):
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

	time_row = "%s\t%s\t%s\t%s\n" % (
		"time", "time", "s", time_string
	)

	# Write line to dynamics file
	dynamics_file.write(time_row)


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
		elif node_id not in duplicate_ids and node_id != "global":
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
		dynamics_string = ", ".join("{0:d}".format(val) for val in dynamics)

	elif datatype == "float":
		dynamics_string = ", ".join("{0:.{1}g}".format(val, DYNAMICS_PRECISION) for val in dynamics)

	elif datatype == "prob":
		dynamics_string = ", ".join("{0:.{1}f}".format(val, PROBABILITY_PRECISION) for val in dynamics)

	elif datatype == "time":
		dynamics_string = ", ".join("{0:.{1}f}".format(val, TIME_PRECISION) for val in dynamics)

	else:
		dynamics_string = dynamics

	return dynamics_string


def read_pathway_file():
	"""
	Reads the pathway file whose filename is specified in PATHWAYS_FILENAME.
	The file is assumed to have 4 columns - pathway IDs, pathway names, the
	list of genes associated with the pathway, and the list of reactions
	associated with the pathway. Note: pathway IDs are currently not being
	read here.
	"""
	pathway_to_genes = {}
	pathway_to_rxns = {}

	with open(PATHWAYS_FILENAME) as pathway_file:
		tsv_reader = csv.reader(pathway_file, delimiter="\t")
		next(tsv_reader, None)  # Ignore header

		# Loop through each row and build dictionary
		for row in tsv_reader:
			pathway_to_genes[row[1]] = eval(row[2])
			pathway_to_rxns[row[1]] = eval(row[3])

	return pathway_to_genes, pathway_to_rxns


def get_pathway_to_nodes(simData, simOutDirs, pathway_to_genes, pathway_to_rxns):
	"""
	Reads simData and constructs dictionary that links each pathway to a set of
	all node IDs that are part of the pathway, starting from the list of
	associated gene and reaction nodes that are given as inputs.
	"""
	# Get bulkMolecule IDs from first simOut directory
	simOutDir = simOutDirs[0]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get reaction IDs from first simOut directory
	simOutDir = simOutDirs[0]
	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = fbaResults.readAttribute("reactionIDs")

	# Get dictionary of genes IDs to RNA IDs
	gene2rna = {}
	for gene_id, rna_id, _ in simData.process.replication.geneData:
		gene2rna[gene_id] = rna_id + "[c]"

	# Get dictionary of RNA IDs to monomer IDs
	rna2monomer = {}
	for monomer_data in simData.process.translation.monomerData:
		rna2monomer[monomer_data[1]] = monomer_data[0]

	# Get TF to RNA data
	rna2tf = {}
	tfToFC = simData.tfToFC
	for tf, transcriptIDdict in tfToFC.items():
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
				if transcriptID + "[c]" in rna2tf:
					rna2tf[transcriptID + "[c]"].append(tfID)
				else:
					rna2tf[transcriptID + "[c]"] = [tfID]

	# Get reaction to metabolite and enzyme data
	reactionStoich = simData.process.metabolism.reactionStoich
	reactionCatalysts = simData.process.metabolism.reactionCatalysts

	# Get list of complex ids in complexation
	complexation_complex_ids = simData.process.complexation.ids_complexes
	equilibrium_complex_ids = simData.process.equilibrium.ids_complexes

	# Initialize dictionary
	pathway_to_nodes = {}

	# Loop through each pathway and its associated genes and reactions
	for pathway_name, gene_list in pathway_to_genes.items():
		rxn_list = pathway_to_rxns[pathway_name]
		node_list = []

		# Loop through all genes associated with the pathway
		for gene_id in gene_list:

			if gene_id in gene2rna:
				node_list.append(gene_id)

				# Get IDs of RNAs and monomers that are produced from the gene
				rna_id = gene2rna[gene_id]
				transcription_id = gene_id + "_TRANSCRIPTION"
				node_list.extend([transcription_id, rna_id])

				if rna_id in rna2monomer:
					monomer_id = rna2monomer[rna_id]
					translation_id = gene_id + "_TRANSLATION"
					node_list.extend([translation_id, monomer_id])

				# Get IDs of TFs and regulation nodes that regulate the gene
				if rna_id in rna2tf:
					tf_ids = rna2tf[rna_id]
					node_list.extend(tf_ids)

					for tf_id in tf_ids:
						regulation_id = tf_id[:-3] + "_" + gene_id + "_REGULATION"
						node_list.append(regulation_id)

		# Loop through all reactions associated with the pathway
		for rxn_id in rxn_list:
			rxn_id_hits = []
			for modeled_rxn_id in reactionIDs:
				if modeled_rxn_id.startswith(rxn_id):
					rxn_id_hits.append(modeled_rxn_id)

			# Check if the reaction is actually being modeled
			for rxn_id in rxn_id_hits:
				node_list.append(rxn_id)

				# Get all metabolites participating in the reaction
				stoich_dict = reactionStoich[rxn_id]
				for metabolite, _ in stoich_dict.items():
					node_list.append(metabolite)

				# Get enzymes that catalyze the reaction
				catalyst_ids = reactionCatalysts.get(rxn_id, [])

				for catalyst in catalyst_ids:
					if catalyst in complexation_complex_ids + equilibrium_complex_ids:
						# Add the catalyst if it is a complex
						node_list.append(catalyst)

						# Add the complexation reaction that forms the complex
						node_list.append(catalyst[:-3] + "_RXN")

		pathway_to_nodes[pathway_name] = list(set(node_list))

	return pathway_to_nodes


def check_nodes_in_pathways(node_ids, pathway_to_nodes):
	"""
	Check if any node in pathway_to_nodes dictionary does not exist in the
	network we constructed.
	"""
	print("Checking pathway-node key...")

	for pathway, node_list in pathway_to_nodes.items():
		for node_id in node_list:
			if node_id not in node_ids:
				print("Node ID %s in pathway %s not found in the list of node IDs." % (node_id, pathway))


def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile=None, metadata=None):
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
	assert ap.n_generation >= N_GENS

	simData = cPickle.load(open(simDataFile))

	# create dict with id: (name, synonyms)
	names_dict = {}
	name_files = [f for f in os.listdir(NAMES_PATHWAY)]
	for file_name in name_files:
		with open(os.path.join(NAMES_PATHWAY, file_name)) as the_file:

			all_data = [line.replace('"', '').replace('\n', '').replace('\r', '').split('\t') for line in the_file.readlines()]
			header = all_data[0]
			data = all_data[1:]

			id_idx = header.index('Object ID')
			synonym_idx = header.index('Synonyms')

			for row in data:
				if row[synonym_idx]:
					synonyms = row[synonym_idx].split(' // ')
					names_dict[row[id_idx]] = (synonyms[0], synonyms[1:])

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

	# Add global nodes to the node list
	add_global_nodes(simData, simOutDirs, node_list)

	# Add state/process-specific nodes and edges to the node list and edge list
	add_replication_and_genes(simData, simOutDirs, node_list, edge_list, names_dict)
	add_transcription_and_transcripts(simData, simOutDirs, node_list, edge_list, names_dict)
	add_translation_and_monomers(simData, simOutDirs, node_list, edge_list, names_dict)
	add_complexation_and_complexes(simData, simOutDirs, node_list, edge_list, names_dict)
	add_metabolism_and_metabolites(simData, simOutDirs, node_list, edge_list, names_dict)
	add_equilibrium(simData, simOutDirs, node_list, edge_list, names_dict)
	add_regulation(simData, simOutDirs, node_list, edge_list, names_dict)

	if GET_PATHWAY_INDEX:
		pathway_to_genes, pathway_to_rxns = read_pathway_file()
		pathway_to_nodes = get_pathway_to_nodes(simData, simOutDirs, pathway_to_genes, pathway_to_rxns)

	# Check for network sanity (optional)
	if CHECK_SANITY:
		print("Performing sanity check on network...")
		node_ids = find_duplicate_nodes(node_list)
		find_runaway_edges(node_ids, edge_list)

		if GET_PATHWAY_INDEX:
			check_nodes_in_pathways(node_ids, pathway_to_nodes)

		print("Sanity check completed.")

	print("Total number of nodes: %d" % (len(node_list)))
	print("Total number of edges: %d" % (len(edge_list)))

	# Open node/edge list files and dynamics file
	nodelist_file = open(os.path.join(plotOutDir, plotOutFileName + "_nodelist.tsv"), 'w')
	edgelist_file = open(os.path.join(plotOutDir, plotOutFileName + "_edgelist.tsv"), 'w')
	dynamics_file = open(os.path.join(plotOutDir, plotOutFileName + "_dynamics.tsv"), 'w')

	if GET_PATHWAY_INDEX:
		pathwaylist_file = open(os.path.join(plotOutDir, plotOutFileName + "_pathwaylist.tsv"), 'w')

	# Write header rows to each of the files
	nodelist_file.write(NODE_LIST_HEADER)
	edgelist_file.write(EDGE_LIST_HEADER)
	dynamics_file.write(DYNAMICS_HEADER)

	if GET_PATHWAY_INDEX:
		pathwaylist_file.write(PATHWAY_LIST_HEADER)

	# Add time and global dynamics data to dynamics file
	add_time_data(simOutDirs, dynamics_file)

	# Write node, edge list and dynamics data tsv files
	for node in node_list:
		node.write_nodelist(nodelist_file)
		node.write_dynamics(dynamics_file)

	for edge in edge_list:
		edge.write_edgelist(edgelist_file)

	# Write pathway data to pathway file
	if GET_PATHWAY_INDEX:
		for pathway_name, node_ids in pathway_to_nodes.items():
			pathwaylist_file.write("%s\t%s\n" % (pathway_name, ", ".join(node_ids)))


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
