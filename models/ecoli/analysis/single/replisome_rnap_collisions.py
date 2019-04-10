"""
Plots the dynamics of the number of collisions between RNAPs and replisomes at
each timestep.

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

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		rnap_data_reader = TableReader(os.path.join(simOutDir, "RnapData"))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		n_aborted_initiations = rnap_data_reader.readColumn("n_aborted_initiations")
		n_collisions = rnap_data_reader.readColumn("n_collisions")

		# Plot
		plt.figure(figsize = (8.5, 11))

		ax = plt.subplot(2, 1, 1)
		ax.plot(time / 60., n_aborted_initiations)
		ax.set_ylabel("Number of aborted transcription initiation events")

		ax = plt.subplot(2, 1, 2)
		ax.plot(time / 60., n_collisions)
		ax.set_ylabel("Number of collisions between replisomes and transcribing RNAPs")
		ax.set_xlabel("Time (min)")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
