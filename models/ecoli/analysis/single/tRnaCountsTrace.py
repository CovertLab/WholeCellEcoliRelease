"""
Plot tRNA counts

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/19/2017
"""

from __future__ import absolute_import

import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get time
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		# Get tRNA IDs and counts
		sim_data = cPickle.load(open(simDataFile, "rb"))
		isTRna = sim_data.process.transcription.rnaData["isTRna"]
		rnaIds = sim_data.process.transcription.rnaData["id"][isTRna]
		charged_trna_ids = sim_data.process.transcription.charged_trna_names

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		mol_indices = {mol: i for i, mol in enumerate(moleculeIds)}

		uncharged_indices = np.array([mol_indices[moleculeId] for moleculeId in rnaIds], np.int)
		charged_indices = np.array([mol_indices[moleculeId] for moleculeId in charged_trna_ids], np.int)

		bulk_counts = bulkMolecules.readColumn("counts")
		rna_counts = bulk_counts[:, uncharged_indices] + bulk_counts[:, charged_indices]

		# Plot
		fig = plt.figure(figsize = (8.5, 11))
		ax = plt.subplot(1, 1, 1)
		ax.plot(time, rna_counts)
		ax.set_xlim([time[0], time[-1]])
		ax.set_xlabel("Time (s)")
		ax.set_ylabel("Counts of tRNAs")
		ax.spines["right"].set_visible(False)
		ax.spines["top"].set_visible(False)
		ax.tick_params(right = "off", top = "off", which = "both", direction = "out")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
