"""
Generates a .json file containing the dynamic locations of molecules bound to
the chromosome.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/2/19
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os
import json

import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath

# Flags to indicate replisome status
NOT_INITIATED = 0
ELONGATING = 1
HAS_TERMINATED = 2

# Flag to indicate RNAP is terminated at this timestep
LAST_TIMESTEP = -1

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Read replichore lengths from sim_data
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		replication_data_reader = TableReader(
			os.path.join(simOutDir, "ReplicationData"))
		rnap_data_reader = TableReader(
			os.path.join(simOutDir, "RnapData"))

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

		# Get list of unique indexes of each replisome
		fork_unique_index_list = np.unique(
			fork_unique_indexes[~np.isnan(fork_unique_indexes)])
		n_unique_replisomes = len(fork_unique_index_list)

		# Parse data such that one row of column corresponds to one unique
		# replisome. The status array is set to the appropriate flag values
		# depending on the status of the corresponding replisome (column)
		# at the given timestep (row).
		fork_coordinates_parsed = np.zeros((n_timesteps, n_unique_replisomes))
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

		# Crop out full columns of NaNs and replace NaNs to zeros for RNAP data
		rnap_isnan = np.isnan(rnap_coordinates)
		n_nan_columns = (rnap_isnan.sum(axis=0) == n_timesteps).sum()

		rnap_status = np.logical_not(rnap_isnan[:, :-n_nan_columns])
		rnap_coordinates = np.nan_to_num(rnap_coordinates[:, :-n_nan_columns])
		rnap_domain_indexes = np.nan_to_num(rnap_domain_indexes[:, :-n_nan_columns])
		rnap_unique_indexes = np.nan_to_num(rnap_unique_indexes[:, :-n_nan_columns])

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

		# Build dictionary of chromosome data
		chromosome_data = {
			"metadata": metadata,
			"time": [round(t, 2) for t in time],
			"right_replichore_len": replichore_lengths[0],
			"left_replichore_len": replichore_lengths[1],
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
				"next_indexes": rnap_next_indexes.tolist(),
				"status": rnap_status.tolist(),
				}
			}

		# Output dictionary to json file
		with open(os.path.join(plotOutDir, plotOutFileName + ".json"), 'w') as f:
			f.write(json.dumps(chromosome_data))


if __name__ == '__main__':
	Plot().cli()
