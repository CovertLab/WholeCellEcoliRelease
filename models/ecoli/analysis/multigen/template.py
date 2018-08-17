"""
Template for multigen analysis plots

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, 'seedOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)

		for sim_dir in ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))

			# Load data
			time = main_reader.readColumn('time')

		plt.figure()

		### Create Plot ###

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
