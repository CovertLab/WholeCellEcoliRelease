"""
Plots a heatmap that marks the times where each enzyme in the metabolic network
is completely absent from the simulation.
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time

		# Load enzyme counts
		enzyme_ids = sim_data.process.metabolism.catalyst_ids
		(enzyme_counts, ) = read_bulk_molecule_counts(simOutDir, (enzyme_ids, ))

		# Highlight enzyme absence
		fig, ax = plt.subplots(figsize=(10, 160))
		ax.imshow(
			(enzyme_counts == 0).T,
			cmap='Reds', interpolation='nearest', aspect='auto'
			)

		ax.set_xticks([0, len(time) - 1])
		ax.set_xticklabels([0, int(time[-1]/60.)])
		ax.set_xlabel('Time (min)')
		ax.set_yticks(np.arange(len(enzyme_ids)))
		ax.set_yticklabels(enzyme_ids)
		ax.set_title(
			'Times when each enzyme is completely absent from the simulation'
			)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
