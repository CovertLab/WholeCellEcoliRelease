"""
Template for cohort analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		cell_paths = self.ap.get_cells()

		# Load data
		## Simple stacking functions for data from all cells
		names = ['ATP[c]']  # Replace with desired list of names
		time = read_stacked_columns(cell_paths, 'Main', 'time')
		(counts,) = read_stacked_bulk_molecules(cell_paths, (names,))

		## Or iterate on each cell if additional processing is needed
		for sim_dir in cell_paths:
			simOutDir = os.path.join(sim_dir, 'simOut')

			# Listeners used
			main_reader = TableReader(os.path.join(simOutDir, 'Main'))

			# Load data
			time = main_reader.readColumn('time')

			(counts,) = read_bulk_molecule_counts(simOutDir, (names,))

		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
