"""
Plot of charged tRNA

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
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
		growth_limits_reader = TableReader(os.path.join(simOutDir, 'GrowthLimits'))

		# Load data
		time = main_reader.readColumn('time')
		fraction_trna_charged = growth_limits_reader.readColumn('fraction_trna_charged')
		net_charged = growth_limits_reader.readColumn('net_charged')

		plt.figure()

		plt.subplot(2,1,1)
		plt.plot(time, fraction_trna_charged)
		plt.ylabel('Fraction tRNA Charged')

		plt.subplot(2,1,2)
		plt.plot(time, net_charged)
		plt.ylabel('Net tRNAs Charged')

		plt.xlabel('Time (s)')

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
