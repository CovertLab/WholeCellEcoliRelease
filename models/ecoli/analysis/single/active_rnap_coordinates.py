"""
Plots the locations of all replisomes and active RNAPs on the chromosome over
time.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2019
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import collections as mc
import numpy as np
from six.moves import cPickle, range

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		replication_data_reader = TableReader(os.path.join(simOutDir, "ReplicationData"))
		RNAP_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		# Convert time unit to minutes
		time_mins = time / 60.

		# Read replisome data
		replisome_coordinates = replication_data_reader.readColumn("fork_coordinates")
		replisome_unique_indexes = replication_data_reader.readColumn("fork_unique_index")

		# Read RNAP data
		RNAP_coordinates = RNAP_data_reader.readColumn("active_rnap_coordinates")
		RNAP_unique_indexes = RNAP_data_reader.readColumn("active_rnap_unique_indexes")

		def parse_coordinates(coordinates, unique_indexes):
			"""
			Parses coordinates array into a list of numpy arrays with two
			columns, where each array stores the time in the first column and
			the coordinates of a unique molecule in the second column.
			"""
			parsed_coordinates = []
			coordinates_dict = {}
			n_rows, n_columns = unique_indexes.shape

			# Loop through unique_indexes array
			for i in range(n_rows):
				for j in range(n_columns):
					if np.isnan(unique_indexes[i, j]):
						continue

					# If unique index is already keyed, update last timestep
					# index and append coordinates value
					if unique_indexes[i, j] in coordinates_dict:
						coordinates_dict[unique_indexes[i, j]]["last_index"] = i
						coordinates_dict[unique_indexes[i, j]]["coordinates"].append(coordinates[i, j])
					# If unique index is first seen, initialize dictionary
					# with first and last timestep indexes, and a list of
					# coordinates
					else:
						coordinates_dict[unique_indexes[i, j]] = {
							"first_index": i,
							"last_index": i,
							"coordinates": [coordinates[i, j]]
							}

			# Loop through values in dictionary
			for v in coordinates_dict.values():
				# Construct numpy array with time column and coordinates column
				# for each molecule
				parsed_coordinates.append(
					np.column_stack((
						time_mins[v["first_index"]:v["last_index"] + 1],
						np.array(v["coordinates"]))))

			return parsed_coordinates

		parsed_replisome_coordinates = parse_coordinates(
			replisome_coordinates, replisome_unique_indexes)
		parsed_RNAP_coordinates = parse_coordinates(
			RNAP_coordinates, RNAP_unique_indexes)

		# Read collision data
		headon_collision_coordinates = RNAP_data_reader.readColumn(
			"headon_collision_coordinates")
		codirectional_collision_coordinates = RNAP_data_reader.readColumn(
			"codirectional_collision_coordinates")

		# Load replichore lengths
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Plot
		plt.figure(figsize = (70, 150))
		ax = plt.subplot2grid((1, 1), (0, 0))

		# Plot a LineCollection of replisome coordinates
		replisome_lc = mc.LineCollection(
			parsed_replisome_coordinates, colors='k', linewidths=2)
		ax.add_collection(replisome_lc)

		# Plot a LineCollection of RNAP coordinates
		RNAP_lc = mc.LineCollection(
			parsed_RNAP_coordinates, colors='#888888', linewidths=0.5)
		ax.add_collection(RNAP_lc)

		# Scatter plot options for collisions
		headon_params = {
			"marker": "x", "markersize": 10, "linewidth": 0,
			"color": "crimson", "label": "head-on"}
		codirectional_params = {
			"marker": "x", "markersize": 10, "linewidth": 0,
			"color": "darkblue", "label": "co-directional"}

		# Plot coordinates of collisions between RNAPs and replisomes. Skip
		# if there are no collisions (no replication forks)
		if headon_collision_coordinates.shape[1] > 0:
			ax.plot(time_mins, headon_collision_coordinates,
				**headon_params)
		if codirectional_collision_coordinates.shape[1] > 0:
			ax.plot(time_mins, codirectional_collision_coordinates,
				**codirectional_params)

		ax.set_xticks([0, time_mins.max()])
		ax.set_yticks([-replichore_lengths[1], 0, replichore_lengths[0]])
		ax.set_yticklabels(['-terC', 'oriC', '+terC'])
		ax.set_ylim([-replichore_lengths[1], replichore_lengths[0]])
		ax.set_ylabel("Position on chromosome (nt)")
		ax.set_xlim([0, time_mins.max()])
		ax.set_xlabel("Time (min)")

		# Add legends for collision markers
		legend_elements = [
			Line2D([0], [0], **headon_params),
			Line2D([0], [0], **codirectional_params)]
		ax.legend(handles=legend_elements)

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata, extension='.pdf')
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
