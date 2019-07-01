"""
Plot RNAse counts

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/14/2015
"""

from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile, "rb"))

		endoRnaseIds = sim_data.process.rna_decay.endoRnaseIds
		exoRnaseIds = sim_data.moleculeGroups.exoRnaseIds
		RNase_IDS = np.concatenate((endoRnaseIds, exoRnaseIds))

		endoRnase_RnaIDs = sim_data.moleculeGroups.endoRnase_RnaIDs
		exoRnase_RnaIDs = sim_data.moleculeGroups.exoRnase_RnaIDs
		RNase_RnaIDS = np.concatenate((endoRnase_RnaIDs, exoRnase_RnaIDs))
		RNase_IDS = np.concatenate((RNase_IDS, RNase_RnaIDS))

		(rnapRnaCounts,) = read_bulk_molecule_counts(simOutDir, (RNase_IDS,))

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		plt.figure(figsize = (8.5, 11))
		plt.rc('xtick', labelsize=7)
		plt.rc('ytick', labelsize=5)

		count = 0
		count_bis = len(RNase_IDS) / 2
		for subplotIdx in xrange(0, len(RNase_IDS)):
			if not subplotIdx % 2:
				rnapRnaCountsIdx = count
				count += 1
			if subplotIdx % 2:
				rnapRnaCountsIdx = count_bis
				count_bis += 1

			ax = plt.subplot(18, 2, 1 + subplotIdx)

			plt.plot(time / 60., rnapRnaCounts[:, rnapRnaCountsIdx])

			if not subplotIdx >= len(RNase_IDS) - 2:
				frame = plt.gca()
				for xlabel_i in frame.axes.get_xticklines():
					xlabel_i.set_visible(True)
				for xlabel_i in frame.axes.get_xticklabels():
					xlabel_i.set_visible(False)

			if subplotIdx >= len(RNase_IDS) - 2:
				plt.xlabel("Time (min)", fontsize = 7)

			if not subplotIdx % 2:
				plt.ylabel("Protein counts", fontsize = 5)
			if subplotIdx % 2:
				plt.ylabel("RNA counts", fontsize = 5)

			plt.title(RNase_IDS[rnapRnaCountsIdx], fontsize = 7)

			max_yticks = 4
			yloc = plt.MaxNLocator(max_yticks)
			ax.yaxis.set_major_locator(yloc)

		plt.subplots_adjust(hspace = 0.75, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

if __name__ == "__main__":
	Plot().cli()
