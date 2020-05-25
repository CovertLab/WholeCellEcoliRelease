"""
Template for cohort analysis plots

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		ap = AnalysisPaths(variantDir, cohort_plot=True)

		for sim_dir in ap.get_cells():
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))

			# Load data
			time = main_reader.readColumn('time')

			names = ['ATP[c]']  # Replace with desired list of names
			(counts,) = read_bulk_molecule_counts(simOutDir, (names,))

		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
