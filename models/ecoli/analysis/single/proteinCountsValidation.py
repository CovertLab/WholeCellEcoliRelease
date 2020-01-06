"""
Compare protein counts to Wisniewski 2014 and Schmidt 2015 data sets

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/3/2015
"""

from __future__ import absolute_import
from __future__ import division

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


def get_sim_wisniewski_counts(
	monomer_counts, wisniewski_monomer_ids, sim_monomer_ids
):
	"""Get simulated monomer counts that correspond to Wisniewski data

	Arguments:
		monomer_counts: Simulated monomer counts (from translation
			process).
		wisniewski_monomer_ids: Monomer IDs from Wisniewski data. IDs
			must appear in same order as in Wisniewski monomer counts.
		sim_monomer_ids: IDs of monomers in the same order as
			monomer_counts.

	Returns:
		The simulated counts of the monomers that appear in the
		Wisniewski data, in the same order as in the Wisniewski data.
	"""
	# type: (np.ndarray, np.ndarray, np.ndarray) -> np.ndarray
	sim_ids_lst = sim_monomer_ids.tolist()
	wisniewski_ids_lst = wisniewski_monomer_ids.tolist()
	wisniewski_idx = [sim_ids_lst.index(x) for x in wisniewski_ids_lst]
	avg_counts = monomer_counts.mean(axis=0)
	return avg_counts[wisniewski_idx]


def get_sim_schmidt_counts(
	monomer_counts, schmidt_monomer_ids, sim_monomer_ids
):
	"""Get simulated monomer counts that correspond to Schmidt data

	Arguments:
		monomer_counts: Simulated monomer counts (from translation
			process).
		schmidt_monomer_ids: Monomer IDs from Schmidt data. IDs
			must appear in same order as in Schmidt monomer counts.
		sim_monomer_ids: IDs of monomers in the same order as
			monomer_counts.

	Returns:
		The simulated counts of the monomers that appear in the
		Schmidt data, in the same order as in the Schmidt data.
	"""
	# type: (np.ndarray, np.ndarray, np.ndarray) -> np.ndarray
	sim_ids_lst = sim_monomer_ids.tolist()
	schmidt_ids_lst = schmidt_monomer_ids.tolist()
	schmidt_idx = [sim_ids_lst.index(x) for x in schmidt_ids_lst]
	avg_counts = monomer_counts.mean(axis=0)
	return avg_counts[schmidt_idx]


class Plot(singleAnalysisPlot.SingleAnalysisPlot):

	def do_plot(
		self, simOutDir, plotOutDir, plotOutFileName, simDataFile,
		validationDataFile, metadata
	):
		if not os.path.isdir(simOutDir):
			raise Exception(
				"simOutDir does not currently exist as a directory"
			)

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile, "rb"))
		validation_data = cPickle.load(open(validationDataFile, "rb"))

		sim_monomer_ids = sim_data.process.translation.monomerData["id"]
		wisniewski_ids = validation_data.protein.wisniewski2014Data["monomerId"]
		schmidt_ids = validation_data.protein.schmidt2015Data["monomerId"]

		monomer_counts_reader = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		monomer_counts = monomer_counts_reader.readColumn("monomerCounts")

		sim_wisniewski_counts = get_sim_wisniewski_counts(
			monomer_counts, wisniewski_ids, sim_monomer_ids)
		sim_schmidt_counts = get_sim_schmidt_counts(
			monomer_counts, schmidt_ids, sim_monomer_ids)

		wisniewski_counts = validation_data.protein.wisniewski2014Data["avgCounts"]
		schmidt_counts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		fig, ax = plt.subplots(2, sharey=True, figsize=(8.5, 11))

		# Wisniewski Counts
		ax[0].scatter(
			np.log10(wisniewski_counts + 1),
			np.log10(sim_wisniewski_counts + 1),
			c='w', edgecolor='k', alpha=.7
		)
		ax[0].set_xlabel("log10(Wisniewski 2014 Counts)")
		ax[0].set_title(
			"Pearson r: %0.2f" %
			pearsonr(
				np.log10(sim_wisniewski_counts + 1),
				np.log10(wisniewski_counts + 1)
			)[0]
		)

		# Schmidt Counts
		ax[1].scatter(
			np.log10(schmidt_counts + 1),
			np.log10(sim_schmidt_counts + 1),
			c='w', edgecolor='k', alpha=.7
		)
		ax[1].set_xlabel("log10(Schmidt 2015 Counts)")
		ax[1].set_title(
			"Pearson r: %0.2f" %
			pearsonr(
				np.log10(sim_schmidt_counts + 1), np.log10(schmidt_counts + 1)
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
