"""
Template for parca analysis plots
"""

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = cPickle.load(f)
		with open(sim_data_file, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validation_data_file, 'rb') as f:
			validation_data = cPickle.load(f)

		plt.figure()

		### Create Plot ###

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
