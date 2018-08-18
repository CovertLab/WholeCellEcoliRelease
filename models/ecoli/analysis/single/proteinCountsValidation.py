"""
Compare protein counts to Wisniewski 2014 data set

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/3/2015
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle
from scipy.stats import pearsonr

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile, "rb"))
		validation_data = cPickle.load(open(validationDataFile, "rb"))

		ids_translation = sim_data.process.translation.monomerData["id"].tolist()
		wisniewski_idx = [ids_translation.index(x) for x in validation_data.protein.wisniewski2014Data["monomerId"].tolist()]
		schmidt_idx = [ids_translation.index(x) for x in validation_data.protein.schmidt2015Data["monomerId"].tolist()]

		monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		avgCounts = monomerCounts.readColumn("monomerCounts").mean(axis=0)
		sim_wisniewski_counts = avgCounts[wisniewski_idx]
		sim_schmidt_counts = avgCounts[schmidt_idx]

		wisniewski_counts = validation_data.protein.wisniewski2014Data["avgCounts"]
		schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		fig, ax = plt.subplots(2, sharey=True, figsize = (8.5, 11))

		# Wisniewski Counts
		ax[0].scatter(np.log10(wisniewski_counts + 1), np.log10(sim_wisniewski_counts + 1), c='w', edgecolor = 'k', alpha=.7)
		ax[0].set_xlabel("log10(Wisniewski 2014 Counts)")
		ax[0].set_title("Pearson r: %0.2f" % pearsonr(np.log10(sim_wisniewski_counts + 1), np.log10(wisniewski_counts + 1))[0])

		# Schmidt Counts
		ax[1].scatter(
			np.log10(schmidt_counts + 1),
			np.log10(sim_schmidt_counts + 1),
			c='w', edgecolor = 'k', alpha=.7)
		ax[1].set_xlabel("log10(Schmidt 2015 Counts)")
		ax[1].set_title("Pearson r: %0.2f" % pearsonr(np.log10(sim_schmidt_counts + 1), np.log10(schmidt_counts + 1))[0])

		plt.ylabel("log10(Simulation Average Counts)")
		# NOTE: This Pearson correlation goes up (at the time of writing) about 0.05 if you only
		# include proteins that you have translational efficiencies for
		plt.xlim(xmin=0)
		plt.ylim(ymin=0)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
