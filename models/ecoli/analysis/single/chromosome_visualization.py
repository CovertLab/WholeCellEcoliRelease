"""
Generates a .json file containing the dynamic locations of molecules bound to
the chromosome.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/2/19
"""

from __future__ import absolute_import, division, print_function

import os
import json

import numpy as np
from six.moves import cPickle, range

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


# Flags to indicate replisome status
NOT_INITIATED = 0
ELONGATING = 1
HAS_TERMINATED = 2

# Flag to indicate RNAP is terminated at this timestep
LAST_TIMESTEP = -1

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Read replichore lengths from sim_data
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Read gene coordinates from sim_data
		gene_start_coordinates = sim_data.process.transcription.rna_data['replication_coordinate']
		gene_direction = sim_data.process.transcription.rna_data['direction']
		gene_direction_rescaled = (2 * (gene_direction - 0.5)).astype(np.int64)
		gene_length = sim_data.process.transcription.rna_data['length'].asNumber(units.nt)
		gene_end_coordinates = gene_start_coordinates + np.multiply(
			gene_direction_rescaled, gene_length)

		# Get common names of genes
		gene_names = []
		for gene_id in sim_data.process.transcription.rna_data['gene_id']:
			if gene_id in sim_data.common_names.genes:
				gene_names.append(sim_data.common_names.genes[gene_id][0])
			else:
				gene_names.append(gene_id)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		replication_data_reader = TableReader(
			os.path.join(simOutDir, "ReplicationData"))
		rnap_data_reader = TableReader(
			os.path.join(simOutDir, "RnapData"))
		rna_synth_prob_reader = TableReader(
			os.path.join(simOutDir, "RnaSynthProb"))

		# Load time data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		n_timesteps = len(time)

		# Load replisome attributes
		fork_coordinates = replication_data_reader.readColumn("fork_coordinates")
		fork_domain_indexes = replication_data_reader.readColumn("fork_domains")
		fork_unique_indexes = replication_data_reader.readColumn("fork_unique_index")

		# Load active RNAP attributes
		rnap_coordinates = rnap_data_reader.readColumn("active_rnap_coordinates")
		rnap_domain_indexes = rnap_data_reader.readColumn("active_rnap_domain_indexes")
		rnap_unique_indexes = rnap_data_reader.readColumn("active_rnap_unique_indexes")
		rnap_n_bound_ribosomes = rnap_data_reader.readColumn("active_rnap_n_bound_ribosomes")

		# Load bound TF attributes
		bound_TF_coordinates = rna_synth_prob_reader.readColumn("bound_TF_coordinates")
		bound_TF_domain_indexes = rna_synth_prob_reader.readColumn("bound_TF_domains")

		# Get list of unique indexes of each replisome
		fork_unique_index_list = np.unique(
			fork_unique_indexes[~np.isnan(fork_unique_indexes)])
		n_unique_replisomes = len(fork_unique_index_list)

		# Parse data such that one row of column corresponds to one unique
		# replisome. The status array is set to the appropriate flag values
		# depending on the status of the corresponding replisome (column)
		# at the given timestep (row).
		fork_coordinates_parsed = np.zeros((n_timesteps, n_unique_replisomes), dtype=np.int64)
		fork_status = np.zeros((n_timesteps, n_unique_replisomes), dtype=np.int64)
		fork_domain_indexes_parsed = np.zeros(n_unique_replisomes, dtype=np.int64)

		for mol_idx, unique_idx in enumerate(fork_unique_index_list):
			time_index, col_index = np.where(fork_unique_indexes == unique_idx)
			fork_coordinates_parsed[time_index, mol_idx] = fork_coordinates[time_index, col_index]
			fork_status[time_index, mol_idx] = ELONGATING

			# Domain indexes are static - just take the first value
			fork_domain_indexes_parsed[mol_idx] = fork_domain_indexes[
				time_index[0], col_index[0]]

		for i in range(fork_status.shape[1]):
			elongating_timesteps = np.where(fork_status[:, i] == ELONGATING)[0]
			fork_status[:elongating_timesteps[0], i] = NOT_INITIATED
			fork_status[(elongating_timesteps[-1] + 1):, i] = HAS_TERMINATED

		# Replace NaNs to zeros for RNAP data
		rnap_isnan = np.isnan(rnap_coordinates)

		rnap_status = np.logical_not(rnap_isnan)
		rnap_coordinates = np.nan_to_num(rnap_coordinates).astype(np.int64)
		rnap_domain_indexes = np.nan_to_num(rnap_domain_indexes).astype(np.int64)
		rnap_unique_indexes = np.nan_to_num(rnap_unique_indexes).astype(np.int64)
		rnap_n_bound_ribosomes = np.nan_to_num(rnap_n_bound_ribosomes).astype(np.int64)

		# Build array for quick indexing into the identical RNAP molecule that
		# shares the same unique ID in the next timestep.
		# Index is set to -1 if the RNAP does not exist at that timestep, or
		# is removed from the simulation at the next timestep.
		rnap_next_indexes = np.empty_like(rnap_coordinates, dtype=np.int64)
		rnap_next_indexes.fill(LAST_TIMESTEP)

		for i in range(rnap_next_indexes.shape[0] - 1):
			for j in range(rnap_next_indexes.shape[1]):
				if not rnap_status[i, j]:
					continue  # RNAP does not exist in this timestep
				next_index = np.where(
					rnap_unique_indexes[i, j] == rnap_unique_indexes[i + 1, :]
					)[0]
				if len(next_index) == 0:
					continue  # RNAP is removed at next timestep
				else:
					rnap_next_indexes[i, j] = next_index[0]

		# Identify all promoter sites that bind to TFs
		mask_nan = np.logical_not(np.isnan(bound_TF_domain_indexes))
		all_TF_domain_indexes = bound_TF_domain_indexes[mask_nan].astype(np.int64)
		all_TF_coordinates = bound_TF_coordinates[mask_nan].astype(np.int64)
		time_index = np.where(mask_nan)[0]

		# Build a boolean array for whether each promoter site is bound to a
		# transcription factor at each timestep
		all_unique_TF_sites, site_indexes = np.unique(
			np.vstack((all_TF_domain_indexes, all_TF_coordinates)).T,
			axis=0, return_inverse=True)
		TF_domain_indexes = all_unique_TF_sites[:, 0]
		TF_coordinates = all_unique_TF_sites[:, 1]

		site_index_to_original_index = {}
		for original_index, site_index in enumerate(site_indexes):
			if site_index in site_index_to_original_index:
				site_index_to_original_index[site_index].append(original_index)
			else:
				site_index_to_original_index[site_index] = [original_index]

		TF_bound = np.zeros(
			(all_unique_TF_sites.shape[0], n_timesteps), dtype=np.bool)

		for site_index, original_indexes in site_index_to_original_index.items():
			TF_bound[site_index, time_index[original_indexes]] = True

		# Build dictionary of chromosome data
		chromosome_data = {
			"metadata": metadata,
			"time": [round(t, 2) for t in time],
			"right_replichore_len": int(replichore_lengths[0]),  # json now fails on np.int64
			"left_replichore_len": int(replichore_lengths[1]),
			"gene_names": gene_names,
			"gene_start_coordinates": gene_start_coordinates.tolist(),
			"gene_end_coordinates": gene_end_coordinates.tolist(),
			"replisomes": {
				"coordinates": fork_coordinates_parsed.tolist(),
				"domain_indexes": fork_domain_indexes_parsed.tolist(),
				"flag_not_initiated": NOT_INITIATED,
				"flag_elongating": ELONGATING,
				"flag_has_terminated": HAS_TERMINATED,
				"status": fork_status.tolist(),
				},
			"active_RNAPs": {
				"coordinates": rnap_coordinates.tolist(),
				"domain_indexes": rnap_domain_indexes.tolist(),
				"flag_last_timestep": LAST_TIMESTEP,
				"n_bound_ribosomes": rnap_n_bound_ribosomes.tolist(),
				"next_indexes": rnap_next_indexes.tolist(),
				"status": rnap_status.tolist(),
				},
			"transcription_factors": {
				"coordinates": TF_coordinates.tolist(),
				"domain_indexes": TF_domain_indexes.tolist(),
				"bound": TF_bound.tolist(),
				},
			}

		# Output dictionary to json file
		with open(os.path.join(plotOutDir, plotOutFileName + ".json"), 'w') as f:
			f.write(json.dumps(chromosome_data))


if __name__ == '__main__':
	Plot().cli()
