"""
Template for variant analysis plots

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/2/18
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		for variant in variants:
			# Load modified variant sim_data
			## Consider calculating variant difference and only loading
			## sim_data once above for better performance.
			with open(ap.get_variant_kb(variant), 'rb') as f:
				variant_sim_data = cPickle.load(f)

			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

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


if __name__ == "__main__":
	Plot().cli()
