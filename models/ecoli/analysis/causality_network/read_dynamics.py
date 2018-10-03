#!/usr/bin/env python
"""
Constructs a causal network of simulation components along with dynamics data
from a multi-generation simulation.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/26/2018
"""
from __future__ import absolute_import
from __future__ import division

import cPickle
import numpy as np
import os

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils import units

from models.ecoli.analysis.causality_network.network_components import Node, Edge, NODELIST_FILENAME, EDGELIST_FILENAME, NODE_LIST_HEADER, EDGE_LIST_HEADER


def add_global_nodes(sim_data, simOutDir, node_list):
	"""
	Add global state nodes to the node list.
	"""
	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get index of full chromosome
	full_chrom_idx = moleculeIDs.index(sim_data.moleculeIds.fullChromosome)

	# Extract dynamics data from each generation
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	counts = bulkMolecules.readColumn("counts")

	cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
	cell_volume = ((1.0/sim_data.constants.cellDensity)*(units.fg*cell_mass)).asNumber(units.L)

	# Append to existing array
	fc_counts_array = counts[:, full_chrom_idx].astype(np.int)
	mass_array = cell_mass
	volume_array = cell_volume

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


def add_replication_and_genes(sim_data, simOutDir, node_list, edge_list, names_dict):
	"""
	Add replication process nodes and gene state nodes with dynamics data to
	the node list, and the edges connected to the replication nodes to the edge
	list. - Heejo
	"""
	dntp_ids = sim_data.moleculeGroups.dNtpIds
	ppi_id = "PPI[c]"
	dnap_ids = ['CPLX0-2361[c]', 'CPLX0-3761[c]', 'CPLX0-3925[c]', 'CPLX0-7910[c]']

	rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
	rnaSynthProb = rnaSynthProbReader.readColumn("rnaSynthProb")

	# Loop through all genes (in the order listed in transcription)
	for i, geneId in enumerate(sim_data.process.transcription.rnaData["geneId"]):
		# Initialize a single gene node
		gene_node = Node("State", "Gene")

		# Add attributes to the node
		# Add common name and synonyms
		if geneId in names_dict:
			attr = {'node_id': geneId,
				'name': names_dict[geneId][0],
				'synonyms': names_dict[geneId][1]
				}
		else:
			attr = {'node_id': geneId,
				'name': geneId
				}

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


def add_transcription_and_transcripts(sim_data, simOutDir, node_list, edge_list, names_dict):
	"""
	Add transcription process nodes and transcript state nodes with dynamics
	data to the node list, and edges connected to the transcription nodes to
	the edge list. - Heejo
	"""
	ntp_ids = sim_data.moleculeGroups.ntpIds
	ppi_id = "PPI[c]"
	rnap_id = "APORNAP-CPLX[c]"

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	counts_array = bulkMolecules.readColumn("counts")

	rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
	n_rna_inits = rnapDataReader.readColumn("rnaInitEvent")

	# Loop through all genes (in the order listed in transcription)
	for i, rnaId in enumerate(sim_data.process.transcription.rnaData["id"]):
		geneId = sim_data.process.transcription.rnaData["geneId"][i]
		isMRna = sim_data.process.transcription.rnaData["isMRna"][i]

		# Initialize a single transcript node
		rna_node = Node("State", "RNA")

		# Add attributes to the node
		# Add common name and synonyms

		# remove compartment tag
		rnaId_no_c = rnaId[:-3]
		if isMRna and geneId in names_dict:
			attr = {'node_id': rnaId,
				'name': names_dict[geneId][0]+'_RNA',
				'synonyms': names_dict[geneId][1]
				}
		elif rnaId_no_c in names_dict:
			attr = {'node_id': rnaId,
				'name': names_dict[rnaId_no_c][0],
				'synonyms': names_dict[rnaId_no_c][1]
				}
		else:
			attr = {'node_id': rnaId,
				'name': rnaId
				}

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
		dynamics = {"transcription initiations": list(n_rna_inits[:, i])}
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


def add_translation_and_monomers(sim_data, simOutDir, node_list, edge_list, names_dict):
	"""
	Add translation process nodes and protein (monomer) state nodes with
	dynamics data to the node list, and edges connected to the translation
	nodes to the edge list. - Heejo
	"""
	# Create nodes for amino acids
	aa_ids = sim_data.moleculeGroups.aaIDs
	gtp_id = "GTP[c]"
	gdp_id = "GDP[c]"
	water_id = "WATER[c]"
	ppi_id = "PPI[c]"

	ribosome_subunit_ids = [sim_data.moleculeIds.s30_fullComplex, sim_data.moleculeIds.s50_fullComplex]
	n_avogadro = sim_data.constants.nAvogadro

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	translationResults = TableReader(os.path.join(simOutDir, "RibosomeData"))
	prob_translation_array = translationResults.readColumn('probTranslationPerTranscript')

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	counts_array = bulkMolecules.readColumn("counts")

	# Extract dynamics data from each generation
	cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
	volume_array = ((1.0/sim_data.constants.cellDensity)*(units.fg*cell_mass)).asNumber(units.L)

	# Loop through all translatable genes
	for idx, data in enumerate(sim_data.process.translation.monomerData):
		monomerId = data[0]
		rnaId = data[1]
		geneId = rnaId.split("_RNA[c]")[0]

		# Initialize a single protein node
		protein_node = Node("State", "Protein")

		# Add attributes to the node
		# TODO: Add common name and synonyms
		monomerId_no_c = monomerId[:-3]
		# Add common name, synonyms, molecular mass
		if monomerId_no_c in names_dict:
			attr = {'node_id': monomerId,
				'name': names_dict[monomerId_no_c][0],
				'synonyms': names_dict[monomerId_no_c][1],
				}
		else:
			attr = {'node_id': monomerId,
				'name': monomerId,
				}
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
		dynamics = {'probability of initiating translation': list(prob_translation_array[:, idx])}
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


def add_complexation_and_complexes(sim_data, simOutDir, node_list, edge_list, names_dict):
	"""
	Add complexation process nodes and complex state nodes with dynamics data
	to the node list, and edges connected to the complexation nodes to the edge
	list. - Eran
	"""

	# TODO (Eran) raw_data is here used to get complexation stoichiometry.
	# This can be saved to sim_data and then retrieved here.
	raw_data = KnowledgeBaseEcoli()

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	n_avogadro = sim_data.constants.nAvogadro

	# List of all complex IDs and reaction IDs
	complex_ids = sim_data.process.complexation.ids_complexes + EQUILIBRIUM_COMPLEXES_IN_COMPLEXATION
	reactionIDs = sim_data.process.complexation.ids_reactions

	complexationResults = TableReader(os.path.join(simOutDir, "ComplexationListener"))
	reactions_array = complexationResults.readColumn('reactionRates')

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	counts_array = bulkMolecules.readColumn("counts")

	# Extract dynamics data from each generation
	cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
	volume_array = ((1.0 / sim_data.constants.cellDensity) * (units.fg * cell_mass)).asNumber(units.L)

	# Get complexation stoichiometry from sim_data
	complexStoich = {}
	# TODO (Eran) save complexationReactions in sim_data, so that raw_data
	# won't be needed
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
		complex_no_c = complex[:-3]
		# Add common name, synonyms, molecular mass
		if complex_no_c in names_dict:
			attr = {'node_id': complex,
				'name': names_dict[complex_no_c][0],
				'synonyms': names_dict[complex_no_c][1],
				'constants': {'mass': 0}
				}
		else:
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


def add_metabolism_and_metabolites(sim_data, simOutDir, node_list, edge_list, names_dict):
	"""
	Add metabolism process nodes and metabolite state nodes with dynamics data
	to the node list, add edges connected to the metabolism nodes to the edge
	list. - Gwanggyu
	Note: forward and reverse reactions are represented as separate nodes.
	"""
	# Get reaction list from first simOut directory
	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = fbaResults.readAttribute("reactionIDs")

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	n_avogadro = sim_data.constants.nAvogadro

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	flux_array = fbaResults.readColumn('reactionFluxes')

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	counts_array = bulkMolecules.readColumn("counts")

	# Extract dynamics data from each generation
	cell_mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass")
	volume_array = ((1.0 / sim_data.constants.cellDensity) * (units.fg * cell_mass)).asNumber(units.L)

	# Get all reaction stoichiometry from sim_data
	reactionStoich = sim_data.process.metabolism.reactionStoich

	# Get reaction to catalyst dict from sim_data
	reactionCatalysts = sim_data.process.metabolism.reactionCatalysts

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
		# Add common name, synonyms, molecular mass
		metabolite_no_c = metabolite[:-3]
		if metabolite_no_c in names_dict:
			attr = {'node_id': metabolite,
				'name': names_dict[metabolite_no_c][0],
				'synonyms': names_dict[metabolite_no_c][1],
				'constants': {'mass': 0}
				}
		else:
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


def add_equilibrium(sim_data, simOutDir, node_list, edge_list, names_dict):
	"""
	Add equilibrium nodes with dynamics data to the node list, and add edges
	connected to the equilibrium nodes to the edge list. - Gwanggyu
	"""
	# Get equilibrium-specific data from sim_data
	equilibriumMoleculeIds = sim_data.process.equilibrium.moleculeNames
	equilibriumRxnIds = sim_data.process.equilibrium.rxnIds
	equilibriumStoichMatrix = sim_data.process.equilibrium.stoichMatrix()
	equilibriumRatesFwd = np.array(sim_data.process.equilibrium.ratesFwd, dtype=np.float32)
	equilibriumRatesRev = np.array(sim_data.process.equilibrium.ratesRev, dtype=np.float32)

	# Get transcription factor-specific data from sim_data
	recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
	tf_ids = sorted(set([x.split("__")[-1] for x in recruitmentColNames if
		x.split("__")[-1] != "alpha"]))
	tfToTfType = sim_data.process.transcription_regulation.tfToTfType

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	equilibriumResults = TableReader(os.path.join(simOutDir, "EquilibriumListener"))
	reactions_array = equilibriumResults.readColumn('reactionRates')

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	counts_array = bulkMolecules.readColumn("counts")

	rnaSynthProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
	pPromoterBoundArray = rnaSynthProb.readColumn("pPromoterBound")

	# Get IDs of complexes that were already added
	complexation_complex_ids = sim_data.process.complexation.ids_complexes

	# Get list of complex IDs in equilibrium
	equilibrium_complex_ids = sim_data.process.equilibrium.ids_complexes

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

	# Get 2CS-specific data from sim_data
	tcsMoleculeIds = sim_data.process.two_component_system.moleculeNames
	tcsRxnIds = sim_data.process.two_component_system.rxnIds
	tcsStoichMatrix = sim_data.process.two_component_system.stoichMatrix()
	tcsRatesFwd = np.array(sim_data.process.two_component_system.ratesFwd, dtype=np.float32)
	tcsRatesRev = np.array(sim_data.process.two_component_system.ratesRev, dtype=np.float32)

	# Initialize list of complex IDs in 2CS (should need instance variable)
	tcs_complex_ids = []

	# Get lists of monomers that were already added
	monomer_ids = []
	for monomerData in sim_data.process.translation.monomerData:
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


def add_regulation(sim_data, simOutDir, node_list, edge_list, names_dict):
	"""
	Add regulation nodes with dynamics data to the node list, and add edges
	connected to the regulation nodes to the edge list. - Gwanggyu
	"""
	# Get regulation-specific data from sim_data
	tfToFC = sim_data.tfToFC

	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	counts_array = bulkMolecules.readColumn("counts")

	# Get IDs of genes and RNAs
	gene_ids = []
	rna_ids = []

	for gene_id, rna_id, _ in sim_data.process.replication.geneData:
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


def add_time_data(simOutDir, dynamics_file):
	"""
	Add time data to the dynamics file. - Gwanggyu
	"""
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

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


def get_pathway_to_nodes(sim_data, simOutDir, pathway_to_genes, pathway_to_rxns):
	"""
	Reads sim_data and constructs dictionary that links each pathway to a set of
	all node IDs that are part of the pathway, starting from the list of
	associated gene and reaction nodes that are given as inputs.
	"""
	# Get bulkMolecule IDs from first simOut directory
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIDs = bulkMolecules.readAttribute("objectNames")

	# Get reaction IDs from first simOut directory
	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	reactionIDs = fbaResults.readAttribute("reactionIDs")

	# Get dictionary of genes IDs to RNA IDs
	gene2rna = {}
	for gene_id, rna_id, _ in sim_data.process.replication.geneData:
		gene2rna[gene_id] = rna_id + "[c]"

	# Get dictionary of RNA IDs to monomer IDs
	rna2monomer = {}
	for monomer_data in sim_data.process.translation.monomerData:
		rna2monomer[monomer_data[1]] = monomer_data[0]

	# Get TF to RNA data
	rna2tf = {}
	tfToFC = sim_data.tfToFC
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
	reactionStoich = sim_data.process.metabolism.reactionStoich
	reactionCatalysts = sim_data.process.metabolism.reactionCatalysts

	# Get list of complex ids in complexation
	complexation_complex_ids = sim_data.process.complexation.ids_complexes
	equilibrium_complex_ids = sim_data.process.equilibrium.ids_complexes

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


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		# Initialize node list and edge list
		node_list = []
		edge_list = []

		# Add global nodes to the node list
		add_global_nodes(sim_data, simOutDir, node_list)

		# Add state/process-specific nodes and edges to the node list and edge list
		add_replication_and_genes(sim_data, simOutDir, node_list, edge_list, names_dict)
		add_transcription_and_transcripts(sim_data, simOutDir, node_list, edge_list, names_dict)
		add_translation_and_monomers(sim_data, simOutDir, node_list, edge_list, names_dict)
		add_complexation_and_complexes(sim_data, simOutDir, node_list, edge_list, names_dict)
		add_metabolism_and_metabolites(sim_data, simOutDir, node_list, edge_list, names_dict)
		add_equilibrium(sim_data, simOutDir, node_list, edge_list, names_dict)
		add_regulation(sim_data, simOutDir, node_list, edge_list, names_dict)

		if GET_PATHWAY_INDEX:
			pathway_to_genes, pathway_to_rxns = read_pathway_file()
			pathway_to_nodes = get_pathway_to_nodes(sim_data, simOutDir, pathway_to_genes, pathway_to_rxns)

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
		dynamics_file = open(os.path.join(plotOutDir, plotOutFileName + "_dynamics.tsv"), 'w')

		# Write header row to dynamics file
		dynamics_file.write(DYNAMICS_HEADER)

		# Add time and global dynamics data to dynamics file
		add_time_data(simOutDir, dynamics_file)

		# Write node, edge list and dynamics data tsv files
		for node in node_list:
			node.write_dynamics(dynamics_file)


if __name__ == '__main__':
	Plot().cli()
