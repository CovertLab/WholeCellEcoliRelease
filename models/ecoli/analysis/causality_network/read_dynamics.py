#!/usr/bin/env python
"""
Reads dynamics data for each of the nodes of a causality network from a single
simulation.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/26/2018
"""
from __future__ import absolute_import, division, print_function

import cPickle
import os
import json
import hashlib

from models.ecoli.analysis import causalityNetworkAnalysis
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils import units

from models.ecoli.analysis.causality_network.network_components import (
	Node, COUNT_UNITS, PROB_UNITS
	)
from models.ecoli.analysis.causality_network.build_network import NODE_ID_SUFFIX
from models.ecoli.processes.metabolism import (
	COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
	)

REQUIRED_COLUMNS = [
	("BulkMolecules", "counts"),
	("ComplexationListener", "complexationEvents"),
	("EquilibriumListener", "reactionRates"),
	("FBAResults", "reactionFluxes"),
	("Mass", "cellMass"),
	("Mass", "dryMass"),
	("Main", "time"),
	("RnaSynthProb", "pPromoterBound"),
	("RnaSynthProb", "rnaSynthProb"),
	("RnaSynthProb", "gene_copy_number"),
	("RnaSynthProb", "n_bound_TF_per_TU"),
	("RnapData", "rnaInitEvent"),
	("RibosomeData", "probTranslationPerTranscript"),
	]

def get_safe_name(s):
	fname = str(int(hashlib.sha256(s.encode('utf-8')).hexdigest(), 16) % 10 **16)
	return fname

class Plot(causalityNetworkAnalysis.CausalityNetworkAnalysis):
	def do_plot(self, simOutDir, seriesOutDir, dynamicsFileName, simDataFile, nodeListFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		filepath.makedirs(seriesOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Read all required tables from simOutDir
		columns = {}
		for table_name, column_name in REQUIRED_COLUMNS:
			columns[(table_name, column_name)] = TableReader(
				os.path.join(simOutDir, table_name)).readColumn(column_name)

		# Convert units of metabolic fluxes in listener to mmol/gCDW/h
		conversion_coeffs = (
			columns[("Mass", "dryMass")] / columns[("Mass", "cellMass")]
			* sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)
			)

		columns[("FBAResults", "reactionFluxesConverted")] = (
			(COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (
			columns[("FBAResults", "reactionFluxes")].T / conversion_coeffs).T
			).asNumber(units.mmol/units.g/units.h)

		# Reshape array for number of bound transcription factors
		n_TU = len(sim_data.process.transcription.rnaData["id"])
		n_TF = len(sim_data.process.transcription_regulation.tf_ids)

		columns[("RnaSynthProb", "n_bound_TF_per_TU")] = (
			columns[("RnaSynthProb", "n_bound_TF_per_TU")]
			).reshape(-1, n_TU, n_TF)

		# Construct dictionaries of indexes where needed
		indexes = {}

		def build_index_dict(id_array):
			return {mol: i for i, mol in enumerate(id_array)}

		molecule_ids = TableReader(
			os.path.join(simOutDir, "BulkMolecules")).readAttribute("objectNames")
		indexes["BulkMolecules"] = build_index_dict(molecule_ids)

		gene_ids = sim_data.process.transcription.rnaData["geneId"]
		indexes["Genes"] = build_index_dict(gene_ids)

		translated_rna_ids = sim_data.process.translation.monomerData["rnaId"]
		indexes["TranslatedRnas"] = build_index_dict(translated_rna_ids)

		metabolism_rxn_ids = TableReader(
			os.path.join(simOutDir, "FBAResults")).readAttribute("reactionIDs")
		indexes["MetabolismReactions"] = build_index_dict(metabolism_rxn_ids)

		complexation_rxn_ids = sim_data.process.complexation.ids_reactions
		indexes["ComplexationReactions"] = build_index_dict(complexation_rxn_ids)

		equilibrium_rxn_ids = sim_data.process.equilibrium.rxnIds
		indexes["EquilibriumReactions"] = build_index_dict(equilibrium_rxn_ids)

		tf_ids = sim_data.process.transcription_regulation.tf_ids
		indexes["TranscriptionFactors"] = build_index_dict(tf_ids)

		# Cache cell volume array (used for calculating concentrations)
		volume = ((1.0 / sim_data.constants.cellDensity) * (
			units.fg * columns[("Mass", "cellMass")])).asNumber(units.L)

		with open(nodeListFile, 'r') as node_file:
			node_dicts = json.load(node_file)

		def dynamics_mapping(dynamics, safe):
			return [{
				'index': index,
				'units': dyn['units'],
				'type': dyn['type'],
				'filename': safe + '.json'}
				for index, dyn in enumerate(dynamics)]

		name_mapping = {}

		def build_dynamics(node_dict):
			node = Node()
			node.node_id = node_dict['ID']
			node.node_type = node_dict['type']
			reader = TYPE_TO_READER_FUNCTION.get(node.node_type)
			if reader:
				reader(sim_data, node, node.node_id, columns, indexes, volume)
			return node

		nodes = [build_dynamics(node_dict) for node_dict in node_dicts]
		nodes.append(time_node(columns))

		for node in nodes:
			dynamics_path = get_safe_name(node.node_id)
			dynamics = node.dynamics_dict()
			dynamics_json = json.dumps(dynamics)

			with open(os.path.join(seriesOutDir, dynamics_path + '.json'), 'w') as dynamics_file:
				dynamics_file.write(dynamics_json)

			name_mapping[node.node_id] = dynamics_mapping(dynamics, dynamics_path)

		root = os.path.dirname(nodeListFile)
		with open(os.path.join(root, 'series.json'), 'w') as series_file:
			series_file.write(json.dumps(name_mapping))


def time_node(columns):
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
	return time_node


def read_global_dynamics(sim_data, node, node_id, columns, indexes, volume):
	"""
	Reads global dynamics from simulation output.
	"""
	cell_mass = columns[("Mass", "cellMass")]

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

	dynamics = {
		"transcription probability": columns[("RnaSynthProb", "rnaSynthProb")][:, gene_index],
		"gene copy number": columns[("RnaSynthProb", "gene_copy_number")][:, gene_index],
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
	count_index = indexes["BulkMolecules"][node_id]

	dynamics = {
		"counts": columns[("BulkMolecules", "counts")][:, count_index],
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
	except Exception:
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
		'flux': columns[("FBAResults", "reactionFluxesConverted")][:, reaction_idx],
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
	except Exception:
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
	tf_idx = indexes["TranscriptionFactors"][tf_id]

	dynamics = {
		'bound TFs': columns[("RnaSynthProb", "n_bound_TF_per_TU")][
			:, gene_idx, tf_idx],
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
