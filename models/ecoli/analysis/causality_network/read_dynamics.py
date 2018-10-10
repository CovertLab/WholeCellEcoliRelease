#!/usr/bin/env python
"""
Reads dynamics data for each of the nodes of a causality network from a single
simulation.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/26/2018
"""
from __future__ import absolute_import, division, print_function

import cPickle
import numpy as np
import os

from models.ecoli.analysis import causalityNetworkAnalysis
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils import units

from models.ecoli.analysis.causality_network.network_components import (
	Node, DYNAMICS_HEADER, COUNT_UNITS, PROB_UNITS
	)
from models.ecoli.analysis.causality_network.build_network import NODE_ID_SUFFIX

REQUIRED_COLUMNS = [
	("BulkMolecules", "counts"),
	("ComplexationListener", "complexationEvents"),
	("EquilibriumListener", "reactionRates"),
	("FBAResults", "reactionFluxes"),
	("Mass", "cellMass"),
	("Main", "time"),
	("RnaSynthProb", "pPromoterBound"),
	("RnaSynthProb", "rnaSynthProb"),
	("RnapData", "rnaInitEvent"),
	("RibosomeData", "probTranslationPerTranscript"),
	]

class Plot(causalityNetworkAnalysis.CausalityNetworkAnalysis):
	def do_plot(self, simOutDir, plotOutDir, dynamicsFileName, simDataFile, nodeListFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Read all required tables from simOutDir
		columns = {}
		for table_name, column_name in REQUIRED_COLUMNS:
			columns[(table_name, column_name)] = TableReader(
				os.path.join(simOutDir, table_name)).readColumn(column_name)

		# Construct dictionaries of indexes where needed
		indexes = {}

		def build_index_dict(id_array):
			return {mol: i for i, mol in enumerate(id_array)}

		molecule_ids = TableReader(
			os.path.join(simOutDir, "BulkMolecules")).readAttribute("objectNames")
		indexes["BulkMolecules"] = build_index_dict(molecule_ids)

		gene_ids = sim_data.process.transcription.rnaData["geneId"]
		indexes["Genes"] = build_index_dict(gene_ids)

		rna_ids = sim_data.process.transcription.rnaData["id"]
		indexes["Rnas"] = build_index_dict(rna_ids)

		translated_rna_ids = sim_data.process.translation.monomerData["rnaId"]
		indexes["TranslatedRnas"] = build_index_dict(translated_rna_ids)

		metabolism_rxn_ids = sim_data.process.metabolism.reactionStoich.keys()
		indexes["MetabolismReactions"] = build_index_dict(metabolism_rxn_ids)

		complexation_rxn_ids = sim_data.process.complexation.ids_reactions
		indexes["ComplexationReactions"] = build_index_dict(complexation_rxn_ids)

		equilibrium_rxn_ids = sim_data.process.equilibrium.rxnIds
		indexes["EquilibriumReactions"] = build_index_dict(equilibrium_rxn_ids)
		
		# Cache cell volume array (used for calculating concentrations)
		volume = ((1.0 / sim_data.constants.cellDensity) * (
			units.fg * columns[("Mass", "cellMass")])).asNumber(units.L)

		# Import attributes of each node from existing node list file
		with open(nodeListFile, 'r') as node_file, open(os.path.join(plotOutDir, dynamicsFileName), 'w') as dynamics_file:
			next(node_file)  # Skip header of node list

			# Add header and time row for dynamics file
			dynamics_file.write(DYNAMICS_HEADER + "\n")
			add_time_row(columns, dynamics_file)

			for line in node_file:
				node = Node()
				node_id, node_type = node.read_attributes_from_tsv(line)

				read_func = TYPE_TO_READER_FUNCTION[node_type]
				read_func(sim_data, node, node_id, columns, indexes, volume)

				node.write_dynamics(dynamics_file)


def add_time_row(columns, dynamics_file):
	"""
	Adds a time row to the dynamics file.
	"""
	time_node = Node()
	attr = {
		'node_class': 'time',
		'node_type': 'time',
		'node_id': 'time',
		}
	time_node.read_attributes(**attr)

	time = columns[("Main", "time")]

	dynamics = {
		'time': time,
		}
	dynamics_units = {
		'time': 's',
		}

	time_node.read_dynamics(dynamics, dynamics_units)
	time_node.write_dynamics(dynamics_file)


def read_global_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads global dynamics from simulation output.
	"""
	cell_mass = columns[("Mass", "cellMass")]

	# Reset IDs of global nodes to "global" in dynamics files
	# (Requested by Fathom)
	node.node_id = "global"

	if node_id == "cell_mass":
		dynamics = {
			"mass": cell_mass,
			}
		dynamics_units = {
			"mass": "fg",
			}

	elif node_id == "cell_volume":
		cell_volume = ((1.0 / sim_data.constants.cellDensity) * (
					units.fg * cell_mass)).asNumber(units.L)
		dynamics = {
			'volume': cell_volume,
			}
		dynamics_units = {
			'volume': 'L',
			}

	else:
		return

	node.read_dynamics(dynamics, dynamics_units)


def read_gene_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for gene nodes from simulation output.
	"""
	gene_index = indexes["Genes"][node_id]
	copy_number_id = sim_data.process.transcription.rnaData["id"][gene_index][:-3] + "__alpha"
	copy_number_index = indexes["BulkMolecules"][copy_number_id]

	dynamics = {
		"transcription probability": columns[("RnaSynthProb", "rnaSynthProb")][:, gene_index],
		"gene copy number": columns[("BulkMolecules", "counts")][:, copy_number_index],
		}
	dynamics_units = {
		"transcription probability": PROB_UNITS,
		"gene copy number": COUNT_UNITS,
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_rna_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for transcript (RNA) nodes from simulation output.
	"""
	rna_index = indexes["Rnas"][node_id]

	dynamics = {
		"counts": columns[("BulkMolecules", "counts")][:, rna_index],
		}
	dynamics_units = {
		"counts": COUNT_UNITS,
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_protein_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for monomer/complex nodes from a simulation output.
	"""
	count_index = indexes["BulkMolecules"][node_id]
	counts = columns[("BulkMolecules", "counts")][:, count_index]
	concentration = (((1 / sim_data.constants.nAvogadro) * counts)/(units.L * volume)).asNumber(units.mmol/units.L)

	dynamics = {
		'counts': counts,
		'concentration': concentration,
		}
	dynamics_units = {
		'counts': COUNT_UNITS,
		'concentration': 'mmol/L',
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_metabolite_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dyanmics data for metabolite nodes from a simulation output.
	"""
	try:
		count_index = indexes["BulkMolecules"][node_id]
	except:
		return  # Metabolite not being modeled
	counts = columns[("BulkMolecules", "counts")][:, count_index]
	concentration = (((1 / sim_data.constants.nAvogadro) * counts)/(units.L * volume)).asNumber(units.mmol/units.L)

	dynamics = {
		'counts': counts,
		'concentration': concentration,
		}
	dynamics_units = {
		'counts': COUNT_UNITS,
		'concentration': 'mmol/L',
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_transcription_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for transcription nodes from simulation output.
	"""
	gene_id = node_id.split(NODE_ID_SUFFIX["transcription"])[0]
	gene_idx = indexes["Genes"][gene_id]

	dynamics = {
		"transcription initiations": columns[("RnapData", "rnaInitEvent")][:, gene_idx],
		}
	dynamics_units = {
		"transcription initiations": COUNT_UNITS,
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_translation_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for translation nodes from a simulation output.
	"""
	rna_id = node_id.split(NODE_ID_SUFFIX["translation"])[0] + "_RNA[c]"
	translation_idx = indexes["TranslatedRnas"][rna_id]

	dynamics = {
		'translation probability': columns[("RibosomeData", "probTranslationPerTranscript")][:, translation_idx],
		}
	dynamics_units = {
		'translation probability': PROB_UNITS,
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_complexation_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for complexation nodes from a simulation output.
	"""
	reaction_idx = indexes["ComplexationReactions"][node_id]

	dynamics = {
		'complexation events': columns[("ComplexationListener", "complexationEvents")][:, reaction_idx],
		}
	dynamics_units = {
		'complexation events': COUNT_UNITS,
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_metabolism_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for metabolism nodes from a simulation output.
	"""
	reaction_idx = indexes["MetabolismReactions"][node_id]

	dynamics = {
		'flux': columns[("FBAResults", "reactionFluxes")][:, reaction_idx],
		}
	dynamics_units = {
		'flux': 'mmol/gCDW/h',
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_equilibrium_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for equilibrium nodes from a simulation output.
	"""
	# TODO (ggsun): Fluxes for 2CS reactions are not being listened to.
	try:
		reaction_idx = indexes["EquilibriumReactions"][node_id]
	except:
		return  # 2CS reaction

	dynamics = {
		'reaction rate': columns[("EquilibriumListener", "reactionRates")][:, reaction_idx],
		}
	dynamics_units = {
		'reaction rate': 'rxns/s',
		}

	node.read_dynamics(dynamics, dynamics_units)


def read_regulation_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads dynamics data for regulation nodes from a simulation output.
	"""
	tf_id, gene_id, _ = node_id.split("_")
	gene_idx = indexes["Genes"][gene_id]
	bound_tf_id = sim_data.process.transcription.rnaData["id"][gene_idx][:-3] + "__" + tf_id
	bound_tf_idx = indexes["BulkMolecules"][bound_tf_id]

	dynamics = {
		'bound TFs': columns[("BulkMolecules", "counts")][:, bound_tf_idx],
		}
	dynamics_units = {
		'bound TFs': COUNT_UNITS,
		}

	node.read_dynamics(dynamics, dynamics_units)


TYPE_TO_READER_FUNCTION = {
	"Global": read_global_dynamics,
	"Gene": read_gene_dynamics,
	"RNA": read_rna_dynamics,
	"Protein": read_protein_dynamics,
	"Complex": read_protein_dynamics,
	"Metabolite": read_metabolite_dynamics,
	"Transcription": read_transcription_dynamics,
	"Translation": read_translation_dynamics,
	"Complexation": read_complexation_dynamics,
	"Equilibrium": read_equilibrium_dynamics,
	"Metabolism": read_metabolism_dynamics,
	"Regulation": read_regulation_dynamics,
	}


if __name__ == '__main__':
	Plot().cli()
