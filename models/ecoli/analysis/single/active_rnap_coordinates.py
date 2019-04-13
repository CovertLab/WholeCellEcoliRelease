"""
Plots the locations of all replisomes and active RNAPs on the chromosome over
time.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2019
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		replication_data_reader = TableReader(os.path.join(simOutDir, "ReplicationData"))
		rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		replisome_coordinates = replication_data_reader.readColumn("fork_coordinates")
		rnap_coordinates = rnap_data_reader.readColumn("active_rnap_coordinates")
		headon_collision_coordinates = rnap_data_reader.readColumn(
			"headon_collision_coordinates")
		codirectional_collision_coordinates = rnap_data_reader.readColumn(
			"codirectional_collision_coordinates")

		# Load replichore lengths
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Plot
		plt.figure(figsize = (20, 150))

		# Scatter plot options for collisions
		headon_params = {
			"marker": "x", "markersize": 10, "linewidth": 0, "color": 'r',
			"label": "head-on"}
		codirectional_params = {
			"marker": "x", "markersize": 10, "linewidth": 0, "color": 'b',
			"label": "co-directional"}

		# Plot coordinates of RNAPs
		plt.plot(time / 60., rnap_coordinates, marker='.', markersize=1,
			linewidth=0, color='gray')

		# Plot coordinates of replisomes
		plt.plot(time / 60., replisome_coordinates, marker='.', markersize=5,
			linewidth=0, color='k')

		# Plot coordinates of collisions between RNAPs and replisomes
		plt.plot(time / 60., headon_collision_coordinates,
			**headon_params)
		plt.plot(time / 60., codirectional_collision_coordinates,
			**codirectional_params)

		plt.xticks([0, time.max() / 60])
		plt.yticks([-replichore_lengths[1], 0, replichore_lengths[0]],
			['-terC', 'oriC', '+terC'])
		plt.ylim([-replichore_lengths[1], replichore_lengths[0]])
		plt.ylabel("Position on chromosome (nt)")
		plt.xlim([0, time.max() / 60])
		plt.xlabel("Time (min)")

		# Add legends for collision markers
		legend_elements = [
			Line2D([0], [0], **headon_params),
			Line2D([0], [0], **codirectional_params)]
		plt.legend(handles=legend_elements)

		plt.tight_layout()

		# Save PNG image (SVG or pdf files are not produced because the file
		# size is too big)
		plt.savefig(os.path.join(plotOutDir, plotOutFileName + '.png'), dpi=200)

		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
