"""
Shows fold change of metabolites over the course of the simulation
"""

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		aa_ids = sim_data.molecule_groups.amino_acids

		# Listeners used
		enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		main_reader = TableReader(os.path.join(simOutDir, "Main"))

		# Metabolite data
		metaboliteNames = np.array(enzymeKineticsdata.readAttribute("metaboliteNames"))
		metaboliteCounts = enzymeKineticsdata.readColumn("metaboliteCountsFinal")
		normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]
		aa_mask = np.array([m in aa_ids for m in metaboliteNames])

		# Highlight outliers
		mean_final = normalizedCounts[-1, ~aa_mask].mean()
		std_final = normalizedCounts[-1, ~aa_mask].std()
		highlighted = (
			(normalizedCounts[-1, :] > mean_final + 3*std_final)
			| (normalizedCounts[-1, :] < mean_final - 3*std_final)
		)

		# Sort amino acids for labeling
		sorted_idx = np.argsort(normalizedCounts[-1, aa_mask])[::-1]

		# Read time info from the listener
		initialTime = main_reader.readAttribute("initialTime")
		time = (main_reader.readColumn("time") - initialTime) / 60

		colors = COLORS_LARGE
		plt.figure(figsize = (8.5, 11))

		# Plot everything but amino acids
		ax = plt.subplot(2, 1, 1)
		ax.set_prop_cycle('color', colors)

		## Plot and label metabolites that are different from the mean
		mask = ~aa_mask & highlighted
		if np.any(mask):
			plt.plot(time, normalizedCounts[:, mask])
			plt.legend(metaboliteNames[mask], fontsize=8)

		## Plot the rest of the metabolites that are not amino acids
		mask = ~aa_mask & ~highlighted
		if np.any(mask):
			plt.plot(time, normalizedCounts[:, mask])

		## Formatting
		plt.xlabel("Time (min)")
		plt.ylabel("Metabolite fold change")
		plt.title('All metabolites (excluding amino acids)')

		# Plot only amino acids
		ax = plt.subplot(2, 1, 2)
		ax.set_prop_cycle('color', colors)
		plt.plot(time, normalizedCounts[:, aa_mask][:, sorted_idx])
		plt.legend(metaboliteNames[aa_mask][sorted_idx], fontsize=8, ncol=2)
		plt.xlabel("Time (min)")
		plt.ylabel("Metabolite fold change - amino acids only")
		plt.title('Only amino acids')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
