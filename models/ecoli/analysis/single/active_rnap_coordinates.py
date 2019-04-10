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

		# Load replichore lengths
		replichore_lengths = sim_data.process.replication.replichore_lengths

		# Plot
		plt.figure(figsize = (20, 150))

		# Plot coordinates of RNAPs
		plt.plot(time / 60., rnap_coordinates, marker='.', markersize=1,
			linewidth=0, color='gray', label="Active RNAPs")

		# Plot coordinates of replisomes
		plt.plot(time / 60., replisome_coordinates, marker='.', markersize=5,
			linewidth=0, color='k', label="Replisomes")

		plt.xticks([0, time.max() / 60])
		plt.yticks([-replichore_lengths[1], 0, replichore_lengths[0]],
			['-terC', 'oriC', '+terC'])
		plt.ylim([-replichore_lengths[1], replichore_lengths[0]])
		plt.ylabel("Position on chromosome (nt)")
		plt.xlim([0, time.max() / 60])
		plt.xlabel("Time (min)")

		plt.tight_layout()

		# Save PNG image (SVG or pdf files are not produced because the file
		# size is too big)
		plt.savefig(os.path.join(plotOutDir, plotOutFileName + '.png'), dpi=200)

		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
