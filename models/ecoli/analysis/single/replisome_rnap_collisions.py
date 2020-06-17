"""
Plots the dynamics of the number of collisions between RNAPs and replisomes at
each timestep.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2019
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		n_total_collisions = rnap_data_reader.readColumn("n_total_collisions")
		n_headon_collisions = rnap_data_reader.readColumn("n_headon_collisions")
		n_codirectional_collisions = rnap_data_reader.readColumn("n_codirectional_collisions")
		n_removed_ribosomes = rnap_data_reader.readColumn("n_removed_ribosomes")

		# Get cumulative sums
		n_total_collisions_cumulative = np.cumsum(n_total_collisions)
		n_headon_collisions_cumulative = np.cumsum(n_headon_collisions)
		n_codirectional_collisions_cumulative = np.cumsum(n_codirectional_collisions)
		n_removed_ribosomes_cumulative = np.cumsum(n_removed_ribosomes)

		# Plot
		plt.figure(figsize = (8.5, 11))

		ax = plt.subplot(4, 1, 1)
		ax.plot(time / 60., n_total_collisions_cumulative)
		ax.set_title("All collisions (Total {})".format(
			n_total_collisions_cumulative[-1]))
		ax.set_ylabel("Cumulative Counts")

		ax = plt.subplot(4, 1, 2)
		ax.plot(time / 60., n_headon_collisions_cumulative)
		ax.set_title("Head-on collisions (Total {})".format(
			n_headon_collisions_cumulative[-1]))
		ax.set_ylabel("Cumulative Counts")

		ax = plt.subplot(4, 1, 3)
		ax.plot(time / 60., n_codirectional_collisions_cumulative)
		ax.set_title("Co-directional collisions (Total {})".format(
			n_codirectional_collisions_cumulative[-1]))
		ax.set_ylabel("Cumulative Counts")

		ax = plt.subplot(4, 1, 4)
		ax.plot(time / 60., n_removed_ribosomes_cumulative)
		ax.set_title("Removed ribosomes (Total {})".format(
			n_removed_ribosomes_cumulative[-1]))
		ax.set_ylabel("Cumulative Counts")
		ax.set_xlabel("Time (min)")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
