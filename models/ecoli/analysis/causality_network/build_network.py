"""
BuildNetwork

Constructs a network representations of simulation components from sim_data,
and generates files for node lists and edge lists.

Adding new nodes to the network:
-------------------------------

To add a new type of nodes to the network (either a state or process), you need
to write a new function within this file (build_network.py), which goes through
all of the instances of the new node type, and for each instance creates a
node:

	new_node = Node()

adds attributes (**attr), which include "node_class", "node_type", "node_id",
"name", and "synonyms":

	new_node.read_attributes(**attr)

and appends the node to the node list:

	self.node_list.append(new_node)

The relevant edges that connect the new node to other nodes also need to be
specified:

	new_edge = Edge("Edge Type")

The source and destination ids for that edge are added with an attribute:

	attr = {
		'src_id': source_id,
		'dst_id': destination_id,
		}

	new_edge.read_attributes(**attr)

and the edge is then added to the edge list:

	self.edge_list.append(new_edge)

With a complete node and edge list, you are ready to add dynamics data to each
node. This is done in read_dynamics.py. You first need to choose appropriate
dynamics data to represents that node's activity, and make sure it is saved in
a listener. read_dynamics.py uses saved listener output to load dynamics into
each node.

read_dynamics.py might require a new function to the dynamics data if it is of
a new node type, specified in the TYPE_TO_READER_FUNCTION dictionary. When the
node list is read, nodes of the new type will be passed into the new function,
which assigns that node dynamics from listener output:

	node.read_dynamics(dynamics, dynamics_units)

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/26/2018

"""
from __future__ import absolute_import, division, print_function

import cPickle
import numpy as np
import os
from itertools import izip

from models.ecoli.analysis.causality_network.network_components import (
	Node, Edge, NODELIST_FILENAME, EDGELIST_FILENAME, NODE_LIST_HEADER,
	EDGE_LIST_HEADER
	)

# Suffixes that are added to the node IDs of a particular type of node
NODE_ID_SUFFIX = {
	"transcription": "_TRS",
	"translation": "_TRL",
	"regulation": "_REG",
	}

# URL template
URL_TEMPLATE = "https://ecocyc.org/ECOLI/substring-search?type=NIL&object={0}&"\
			   "quickSearch=Quick+Search"

"""
The following groups of molecules participate in multiple processes and are
thus identified here to prevent the addition of duplicate nodes.

Note:
	Identifying multi-process participatory molecules in this way is not
	required because --check_sanity checks for duplicate nodes. However,
	identifying such molecules here can streamline network building by
	eliminating the need to search through nodes that were added previously.

TODO:
	Future proof by programmatically finding multi-process participatory
	molecules (maybe dictionary of booleans).
"""
# Proteins that are reactants or products of a metabolic reaction
PROTEINS_IN_METABOLISM = ["EG50003-MONOMER[c]", "PHOB-MONOMER[c]",
	"PTSI-MONOMER[c]", "PTSH-MONOMER[c]"]

# Equilibrium complexes that are formed from deactivated equilibrium reactions,
# but are reactants in a complexation reaction
EQUILIBRIUM_COMPLEXES_IN_COMPLEXATION = ["CPLX0-7620[c]", "CPLX0-7701[c]",
	"CPLX0-7677[c]", "MONOMER0-1781[c]", "CPLX0-7702[c]"]

# Metabolites that are used as ligands in equilibrium, but do not participate
# in any metabolic reactions
METABOLITES_ONLY_IN_EQUILIBRIUM = ["4FE-4S[c]", "NITRATE[p]"]

# Molecules in 2CS (two component system) reactions that are not proteins
NONPROTEIN_MOLECULES_IN_2CS = ["ATP[c]", "ADP[c]", "WATER[c]", "PI[c]",
	"PROTON[c]", "PHOSPHO-PHOB[c]"]


class BuildNetwork(object):
	"""
	Constructs a causality network of simulation components, namely states and
	processes, of a whole-cell simulation using sim_data. Writes two files
	(node list and edge list) that are subsequently used by the dynamics reader
	to extract simulation results, and for the visual representation of the
	network.
	"""
	def __init__(self, sim_data_file, output_dir, check_sanity=False):
		"""
		Args:
			sim_data_file: path to the variant sim_data cPickle file used for
			building the network.
			output_dir: output directory for the node list and edge list files.
			check_sanity: if set to True, checks if there are any nodes with
			duplicate IDs in the network.
			TODO: have check_sanity looks for disconnected nodes, and edges
			with non-existent nodes.
		"""
		# Open simulation data and save as attribute
		with open(sim_data_file, 'rb') as f:
			self.sim_data = cPickle.load(f)

		self.output_dir = output_dir
		self.check_sanity = check_sanity

		self.names_dict = self.sim_data.fathom.common_names

		self.node_list = []
		self.edge_list = []


	def run(self):
		"""
		Build the network and write node/edge list files.
		"""
		self._build_network()
		self._write_files()


	def _build_network(self):
		"""
		Add nodes and edges to the node/edge lists, and check for network
		sanity (optional).
		"""
		# Add global nodes to the node list
		self._add_global_nodes()

		# Add state/process-specific nodes and edges to the node and edge list
		self._add_genes()
		self._add_transcription_and_transcripts()
		self._add_translation_and_monomers()
		self._add_complexation_and_complexes()
		self._add_metabolism_and_metabolites()
		self._add_equilibrium()
		self._add_regulation()

		# Check for network sanity (optional)
		if self.check_sanity:
			self._find_duplicate_nodes()


	def _write_files(self):
		"""
		Write node/edge list as separate .tsv files.
		"""
		# Open node list file
		with open(os.path.join(
				self.output_dir, NODELIST_FILENAME), 'w') as nodelist_file:

			# Write header row
			nodelist_file.write(NODE_LIST_HEADER + "\n")

			# Write one row for each node
			for node in self.node_list:
				node.write_nodelist(nodelist_file)

		# Open edge list file
		with open(os.path.join(
				self.output_dir, EDGELIST_FILENAME), 'w') as edgelist_file:

			# Write header row
			edgelist_file.write(EDGE_LIST_HEADER + "\n")

			# Write one row for each edge
			for edge in self.edge_list:
				edge.write_edgelist(edgelist_file)


	def _add_global_nodes(self):
		"""
		Add global state nodes to the node list.
		"""
		# Add total cell mass node to node list
		mass_node = Node()
		attr = {
			"node_class": "State",
			"node_type": "Global",
			"node_id": "cell_mass",
			"name": "Total cell mass"
			}
		mass_node.read_attributes(**attr)

		# Add total cell volume node to node list
		volume_node = Node()
		attr = {
			"node_class": "State",
			"node_type": "Global",
			"node_id": "cell_volume",
			"name": "Total cell volume"
			}
		volume_node.read_attributes(**attr)

		self.node_list.extend([mass_node, volume_node])


	def _add_genes(self):
		"""
		Add gene state nodes to the node list.
		"""
		# Loop through all genes (in the order listed in transcription)
		for gene_id in self.sim_data.process.transcription.rnaData["geneId"]:
			
			# Initialize a single gene node
			gene_node = Node()

			# Get name and synonyms for gene
			gene_name, gene_synonym = self.names_dict.get(gene_id, (gene_id, [gene_id]))

			# Get URL for gene
			gene_url = URL_TEMPLATE.format(gene_id)

			attr = {
				"node_class": "State",
				"node_type": "Gene",
				"node_id": gene_id,
				"name": gene_name,
				"synonyms": gene_synonym,
				"url": gene_url,
				}

			gene_node.read_attributes(**attr)

			# Append gene node to node_list
			self.node_list.append(gene_node)


	def _add_transcription_and_transcripts(self):
		"""
		Add transcription process nodes and transcript state nodes to the node
		list, and edges connected to the transcription nodes to the edge list.
		"""
		ntp_ids = self.sim_data.moleculeGroups.ntpIds
		ppi_id = "PPI[c]"
		rnap_id = self.sim_data.moleculeIds.rnapFull

		# Loop through all genes (in the order listed in transcription)
		for rna_id, gene_id, is_mrna in izip(
				self.sim_data.process.transcription.rnaData["id"],
				self.sim_data.process.transcription.rnaData["geneId"],
				self.sim_data.process.transcription.rnaData["isMRna"]):

			# Initialize a single transcript node
			rna_node = Node()

			# Add attributes to the node
			# Add common name and synonyms
			# Remove compartment tag
			rna_id_no_compartment = rna_id[:-3]
			gene_name, gene_synonyms = self.names_dict.get(gene_id,
				(gene_id, [gene_id]))

			if is_mrna:
				rna_name = gene_name + " mRNA"
				rna_synonyms = [x + " mRNA" for x in gene_synonyms]
			else:
				rna_name, rna_synonyms = self.names_dict.get(rna_id_no_compartment, (rna_id, [rna_id]))

			attr = {
				'node_class': 'State',
				'node_type': 'RNA',
				'node_id': rna_id,
				'name': rna_name,
				'synonyms': rna_synonyms,
				}

			rna_node.read_attributes(**attr)

			# Append transcript node to node_list
			self.node_list.append(rna_node)

			# Initialize a single transcription node for each gene
			transcription_node = Node()

			# Add attributes to the node
			transcription_id = gene_id + NODE_ID_SUFFIX["transcription"]
			transcription_name = gene_name + " transcription"
			transcription_synonyms = [x + " transcription" for x in gene_synonyms]
			attr = {
				'node_class': 'Process',
				'node_type': 'Transcription',
				'node_id': transcription_id,
				'name': transcription_name,
				'synonyms': transcription_synonyms,
				}
			transcription_node.read_attributes(**attr)

			# Append transcription node to node_list
			self.node_list.append(transcription_node)

			# Add edge from gene to transcription node
			self._append_edge("Transcription", gene_id, transcription_id)

			# Add edge from transcription to transcript node
			self._append_edge("Transcription", transcription_id, rna_id)

			# Add edges from NTPs to transcription nodes
			for ntp_id in ntp_ids:
				self._append_edge("Transcription", ntp_id, transcription_id)

			# Add edge from transcription to Ppi
			self._append_edge("Transcription", transcription_id, ppi_id)

			# Add edges from RNA polymerases to transcription
			self._append_edge("Transcription", rnap_id, transcription_id)


	def _add_translation_and_monomers(self):
		"""
		Add translation process nodes and protein (monomer) state nodes to the
		node list, and edges connected to the translation nodes to the edge
		list.
		"""
		# Create nodes for amino acids
		aa_ids = self.sim_data.moleculeGroups.aaIDs
		gtp_id = "GTP[c]"
		gdp_id = "GDP[c]"
		water_id = "WATER[c]"
		ppi_id = "PPI[c]"

		ribosome_subunit_ids = [self.sim_data.moleculeIds.s30_fullComplex,
			self.sim_data.moleculeIds.s50_fullComplex]

		# Construct dictionary to get corrensponding gene IDs from RNA IDs
		rna_id_to_gene_id = {}
		for rna_id, gene_id in izip(
				self.sim_data.process.transcription.rnaData["id"],
				self.sim_data.process.transcription.rnaData["geneId"]):
			rna_id_to_gene_id[rna_id] = gene_id

		# Loop through all translatable genes
		for monomer_id, rna_id in izip(
				self.sim_data.process.translation.monomerData["id"],
				self.sim_data.process.translation.monomerData["rnaId"]):

			gene_id = rna_id_to_gene_id[rna_id]

			# Initialize a single protein node
			protein_node = Node()

			# Add attributes to the node
			monomer_id_no_compartment = monomer_id[:-3]
			monomer_name, monomer_synonyms = self.names_dict.get(monomer_id_no_compartment,
				(monomer_id, [monomer_id]))
			gene_name, gene_synonyms = self.names_dict.get(gene_id,
				(gene_id, [gene_id]))

			attr = {
				'node_class': 'State',
				'node_type': 'Protein',
				'node_id': monomer_id,
				'name': monomer_name,
				'synonyms': monomer_synonyms,
				}

			protein_node.read_attributes(**attr)

			# Append protein node to node_list
			self.node_list.append(protein_node)

			# Initialize a single translation node for each transcript
			translation_node = Node()

			# Add attributes to the node
			translation_id = gene_id + NODE_ID_SUFFIX["translation"]
			translation_name = gene_name + " translation"
			translation_synonyms = [x + " translation" for x in gene_synonyms]
			attr = {
				'node_class': 'Process',
				'node_type': 'Translation',
				'node_id': translation_id,
				'name': translation_name,
				'synonyms': translation_synonyms,
				}
			translation_node.read_attributes(**attr)

			# Append translation node to node_list
			self.node_list.append(translation_node)

			# Add edge from transcript to translation node
			self._append_edge("Translation", rna_id, translation_id)

			# Add edge from translation to monomer node
			self._append_edge("Translation", translation_id, monomer_id)

			# Add edges from amino acids to translation node
			for aa_id in aa_ids:
				self._append_edge("Translation", aa_id, translation_id)

			# Add edges from other reactants to translation node
			for reactant_id in [gtp_id, water_id]:
				self._append_edge("Translation", reactant_id, translation_id)

			# Add edges from translation to other product nodes
			for product_id in [gdp_id, ppi_id, water_id]:
				self._append_edge("Translation", translation_id, product_id)

			# Add edges from ribosome subunits to translation node
			for subunit_id in ribosome_subunit_ids:
				self._append_edge("Translation", subunit_id, translation_id)


	def _add_complexation_and_complexes(self):
		"""
		Add complexation process nodes and complex state nodes to the node
		list, and edges connected to the complexation nodes to the edge
		list.
		"""
		# List of all complex IDs and reaction IDs
		complex_ids = self.sim_data.process.complexation.ids_complexes + EQUILIBRIUM_COMPLEXES_IN_COMPLEXATION
		reaction_ids = self.sim_data.process.complexation.ids_reactions

		molecule_ids = self.sim_data.process.complexation.moleculeNames
		stoich_matrix = self.sim_data.process.complexation.stoichMatrix()

		# Loop through all complexation reactions
		for reaction_index, reaction_id in enumerate(reaction_ids):
			# Initialize a single complexation node for each complexation reaction
			complexation_node = Node()

			# Add attributes to the node
			attr = {
				'node_class': 'Process',
				'node_type': 'Complexation',
				'node_id': reaction_id,
				'name': reaction_id,
				}
			complexation_node.read_attributes(**attr)

			# Append node to node_list
			self.node_list.append(complexation_node)

			# Get reaction stoichiometry from stoichimetric matrix
			stoich_vector = stoich_matrix[:, reaction_index]
			molecule_indices = np.where(stoich_vector)[0]
			stoich_coeffs = stoich_vector[molecule_indices]

			# Loop through all proteins participating in the reaction
			for molecule_index, stoich in izip(molecule_indices, stoich_coeffs):
				# Add complexation edges
				# Note: the direction of the edge is determined by the sign of the
				# stoichiometric coefficient.
				if stoich > 0:
					self._append_edge("Complexation", reaction_id,
						molecule_ids[molecule_index], stoich)
				else:
					self._append_edge("Complexation",
						molecule_ids[molecule_index], reaction_id, stoich)


		for complex_id in complex_ids:
			# Initialize a single complex node for each complex
			complex_node = Node()

			# Add attributes to the node
			complex_id_no_compartment = complex_id[:-3]
			complex_name, complex_synonyms = self.names_dict.get(
				complex_id_no_compartment, (complex_id, [complex_id]))

			attr = {
				'node_class': 'State',
				'node_type': 'Complex',
				'node_id': complex_id,
				'name': complex_name,
				'synonyms': complex_synonyms,
				}

			complex_node.read_attributes(**attr)

			# Append node to node_list
			self.node_list.append(complex_node)


	def _add_metabolism_and_metabolites(self):
		"""
		Add metabolism process nodes and metabolite state nodes to the node
		list, add edges connected to the metabolism nodes to the edge list.
		Note: forward and reverse reactions are represented as separate nodes.
		"""
		# Get all reaction stoichiometry from sim_data
		reaction_stoich = self.sim_data.process.metabolism.reactionStoich

		# Get reaction to catalyst dict from sim_data
		reaction_catalysts = self.sim_data.process.metabolism.reactionCatalysts

		# Initialize list of metabolite IDs
		metabolite_ids = []

		# Loop through all reactions
		for reaction_id, stoich_dict in reaction_stoich.iteritems():
			# Initialize a single metabolism node for each reaction
			metabolism_node = Node()

			# Add attributes to the node
			attr = {
				'node_class': 'Process',
				'node_type': 'Metabolism',
				'node_id': reaction_id,
				'name': reaction_id,
				}
			metabolism_node.read_attributes(**attr)

			# Append node to node_list
			self.node_list.append(metabolism_node)

			# Get list of proteins that catalyze this reaction
			catalyst_list = reaction_catalysts.get(reaction_id, [])

			# Add an edge from each catalyst to the metabolism node
			for catalyst in catalyst_list:
				self._append_edge("Metabolism", catalyst, reaction_id)

			# Loop through all metabolites participating in the reaction
			for metabolite, stoich in stoich_dict.items():

				# Add metabolites that were not encountered
				if metabolite not in metabolite_ids:
					metabolite_ids.append(metabolite)

				# Add Metabolism edges
				# Note: the direction of the edge is determined by the sign of the
				# stoichiometric coefficient.
				if stoich > 0:
					self._append_edge("Metabolism", reaction_id, metabolite,
						stoich)
				else:
					self._append_edge("Metabolism", metabolite, reaction_id,
						stoich)

		# Loop through all metabolites
		for metabolite_id in metabolite_ids:
			# Skip proteins - they should have already been added
			if metabolite_id in PROTEINS_IN_METABOLISM:
				continue

			# Initialize a single metabolite node for each metabolite
			metabolite_node = Node()

			# Add attributes to the node
			metabolite_id_no_compartment = metabolite_id[:-3]
			metabolite_name, metabolite_synonyms = self.names_dict.get(
				metabolite_id_no_compartment, (metabolite_id, [metabolite_id]))

			attr = {
				'node_class': 'State',
				'node_type': 'Metabolite',
				'node_id': metabolite_id,
				'name': metabolite_name,
				'synonyms': metabolite_synonyms,
				}

			metabolite_node.read_attributes(**attr)

			# Append node to node_list
			self.node_list.append(metabolite_node)


	def _add_equilibrium(self):
		"""
		Add equilibrium nodes to the node list, and add edges connected to the
		equilibrium nodes to the edge list.
		"""
		# Get equilibrium-specific data from sim_data
		equilibrium_molecule_ids = self.sim_data.process.equilibrium.moleculeNames
		equilibrium_reaction_ids = self.sim_data.process.equilibrium.rxnIds
		equilibrium_stoich_matrix = self.sim_data.process.equilibrium.stoichMatrix()

		# Get IDs of complexes that were already added
		complexation_complex_ids = self.sim_data.process.complexation.ids_complexes

		# Get list of complex IDs in equilibrium
		equilibrium_complex_ids = self.sim_data.process.equilibrium.ids_complexes

		# Loop through each equilibrium reaction
		for reaction_index, reaction_id in enumerate(equilibrium_reaction_ids):

			# Initialize a single equilibrium node for each equilibrium reaction
			equilibrium_node = Node()

			# Add attributes to the node
			reaction_name = reaction_id[:-4] + " equilibrium rxn"
			attr = {
				'node_class': 'Process',
				'node_type': 'Equilibrium',
				'node_id': reaction_id,
				'name': reaction_name,
				}
			equilibrium_node.read_attributes(**attr)

			# Append new node to node_list
			self.node_list.append(equilibrium_node)

			# Extract column corresponding to reaction in the stoichiometric matrix
			equilibrium_stoich_matrix_column = equilibrium_stoich_matrix[:, reaction_index]

			# Loop through each element in column
			for molecule_index, stoich in enumerate(equilibrium_stoich_matrix_column):
				molecule_id = equilibrium_molecule_ids[molecule_index]

				# Add Equilibrium edges
				# Note: the direction of the edge is determined by the sign of the
				# stoichiometric coefficient.
				if stoich > 0:
					self._append_edge("Equilibrium", reaction_id, molecule_id,
						stoich)
				else:
					self._append_edge("Equilibrium", molecule_id, reaction_id,
						stoich)

		# Get 2CS-specific data from sim_data
		tcs_molecule_ids = self.sim_data.process.two_component_system.moleculeNames
		tcs_reaction_ids = self.sim_data.process.two_component_system.rxnIds
		tcs_stoich_matrix = self.sim_data.process.two_component_system.stoichMatrix()

		# Initialize list of complex IDs in 2CS
		# TODO (ggsun): add this to sim_data
		tcs_complex_ids = []

		# Get lists of monomers that were already added
		monomer_ids = list(self.sim_data.process.translation.monomerData["id"])

		# Loop through each 2CS reaction
		for reaction_index, reaction_id in enumerate(tcs_reaction_ids):

			# Initialize a single equilibrium node for each equilibrium reaction
			equilibrium_node = Node()

			# Add attributes to the node
			reaction_name = reaction_id[:-4] + " 2CS rxn"
			attr = {
				'node_class': 'Process',
				'node_type': 'Equilibrium',
				'node_id': reaction_id,
				'name': reaction_name,
				}
			equilibrium_node.read_attributes(**attr)

			# Append new node to node_list
			self.node_list.append(equilibrium_node)

			# Extract column corresponding to reaction in the stoichiometric matrix
			tcs_stoich_matrix_column = tcs_stoich_matrix[:, reaction_index]

			# Loop through each element in column
			for molecule_index, stoich in enumerate(tcs_stoich_matrix_column):
				molecule_id = tcs_molecule_ids[molecule_index]

				if molecule_id not in monomer_ids + NONPROTEIN_MOLECULES_IN_2CS:
					tcs_complex_ids.append(molecule_id)

				# Add Equilibrium edges
				# Note: the direction of the edge is determined by the sign of the
				# stoichiometric coefficient.
				if stoich > 0:
					self._append_edge("Equilibrium", reaction_id, molecule_id,
						stoich)
				else:
					self._append_edge("Equilibrium", molecule_id, reaction_id,
						stoich)

		# Add new complexes that were encountered here
		for complex_id in list(set(equilibrium_complex_ids + tcs_complex_ids)):
			if complex_id in complexation_complex_ids:
				continue

			# Initialize a single complex node for each complex
			complex_node = Node()

			# Add attributes to the node
			complex_id_no_compartment = complex_id[:-3]
			complex_name, complex_synonyms = self.names_dict.get(
				complex_id_no_compartment, (complex_id, [complex_id]))

			attr = {
				'node_class': 'State',
				'node_type': 'Complex',
				'node_id': complex_id,
				'name': complex_name,
				'synonyms': complex_synonyms,
				}
			complex_node.read_attributes(**attr)

			# Append node to node_list
			self.node_list.append(complex_node)

		# Loop through metabolites that only appear in equilibrium
		for metabolite_id in METABOLITES_ONLY_IN_EQUILIBRIUM:
			# Initialize a single metabolite node for each metabolite
			metabolite_node = Node()

			# Add attributes to the node
			metabolite_id_no_compartment = metabolite_id[:-3]
			metabolite_name, metabolite_synonyms = self.names_dict.get(
				metabolite_id_no_compartment, (metabolite_id, [metabolite_id]))

			attr = {
				'node_class': 'State',
				'node_type': 'Metabolite',
				'node_id': metabolite_id,
				'name': metabolite_name,
				'synonyms': metabolite_synonyms,
				}
			metabolite_node.read_attributes(**attr)

			# Append node to node_list
			self.node_list.append(metabolite_node)


	def _add_regulation(self):
		"""
		Add regulation nodes with to the node list, and add edges connected to
		the regulation nodes to the edge list.
		"""
		# Get recruitment column names from sim_data
		recruitment_column_names = self.sim_data.process.transcription_regulation.recruitmentColNames

		# Get IDs of genes and RNAs
		gene_ids = self.sim_data.process.replication.geneData["name"]
		rna_ids = list(self.sim_data.process.replication.geneData["rnaId"])

		# Loop through all recruitment column names
		for recruitment_column_name in recruitment_column_names:
			if recruitment_column_name.endswith("__alpha"):
				continue

			rna_id, tf = recruitment_column_name.split("__")

			# Add localization ID to the TF ID
			tf_id = tf + "[c]"

			# Find corresponding ID of gene
			gene_id = gene_ids[rna_ids.index(rna_id)]

			# Initialize a single regulation node for each TF-gene pair
			regulation_node = Node()

			# Add attributes to the node
			reg_id = tf + "_" + gene_id + NODE_ID_SUFFIX["regulation"]
			reg_name = tf + "-" + gene_id + " gene regulation"
			attr = {
				'node_class': 'Process',
				'node_type': 'Regulation',
				'node_id': reg_id,
				'name': reg_name,
				}
			regulation_node.read_attributes(**attr)

			self.node_list.append(regulation_node)

			# Add edge from TF to this regulation node
			self._append_edge("Regulation", tf_id, reg_id)

			# Add edge from this regulation node to the gene
			self._append_edge("Regulation", reg_id, gene_id)


	def _find_duplicate_nodes(self):
		"""
		Identify nodes that have duplicate IDs.
		"""
		node_ids = []

		# Loop through all nodes in the node_list
		for node in self.node_list:
			# Get ID of the node
			node_ids.append(node.get_node_id())

		duplicate_ids = set([x for x in node_ids if node_ids.count(x) > 1])

		# Print duplicate node IDs that were found
		if len(duplicate_ids) > 0:
			raise Exception("%d node IDs were found to be duplicate: %s"
				% (len(duplicate_ids), duplicate_ids))


	def _append_edge(self, type, src, dst, stoichiometry=""):
		"""
		Helper function for appending new nodes to the network.
		"""
		edge = Edge(type)
		attr = {
			'src_id': src,
			'dst_id': dst,
			'stoichiometry': stoichiometry,
			}
		edge.read_attributes(**attr)
		self.edge_list.append(edge)
