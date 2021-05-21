"""
Compare protein counts to Wisniewski 2014 and Schmidt 2015 data sets
"""

from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

from wholecell.io.tablereader import TableReader
from wholecell.utils.protein_counts import get_simulated_validation_counts
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):

	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
			validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile, "rb"))
		validation_data = cPickle.load(open(validationDataFile, "rb"))

		sim_monomer_ids = sim_data.process.translation.monomer_data["id"]
		wisniewski_ids = validation_data.protein.wisniewski2014Data["monomerId"]
		schmidt_ids = validation_data.protein.schmidt2015Data["monomerId"]
		wisniewski_counts = validation_data.protein.wisniewski2014Data["avgCounts"]
		schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		monomer_counts = monomer_counts_reader.readColumn("monomerCounts")

		sim_wisniewski_counts, val_wisniewski_counts = get_simulated_validation_counts(
			wisniewski_counts, monomer_counts, wisniewski_ids, sim_monomer_ids)
		sim_schmidt_counts, val_schmidt_counts = get_simulated_validation_counts(
			schmidt_counts, monomer_counts, schmidt_ids, sim_monomer_ids)

		# noinspection PyTypeChecker
		fig, ax = plt.subplots(2, sharey=True, figsize=(8.5, 11))

		# Wisniewski Counts
		ax[0].scatter(
			np.log10(val_wisniewski_counts + 1),
			np.log10(sim_wisniewski_counts + 1),
			c='w', edgecolor='k', alpha=.7
		)
		ax[0].set_xlabel("log10(Wisniewski 2014 Counts)")
		ax[0].set_title(
			"Pearson r: %0.2f" %
			pearsonr(
				np.log10(sim_wisniewski_counts + 1),
				np.log10(val_wisniewski_counts + 1)
			)[0]
		)

		# Schmidt Counts
		ax[1].scatter(
			np.log10(val_schmidt_counts + 1),
			np.log10(sim_schmidt_counts + 1),
			c='w', edgecolor='k', alpha=.7
		)
		ax[1].set_xlabel("log10(Schmidt 2015 Counts)")
		ax[1].set_title(
			"Pearson r: %0.2f" %
			pearsonr(
				np.log10(sim_schmidt_counts + 1), np.log10(val_schmidt_counts + 1)
			)[0]
		)

		plt.ylabel("log10(Simulation Average Counts)")
		# NOTE: This Pearson correlation goes up (at the time of
		# writing) about 0.05 if you only include proteins that you have
		# translational efficiencies for
		plt.xlim(xmin=0)
		plt.ylim(ymin=0)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
