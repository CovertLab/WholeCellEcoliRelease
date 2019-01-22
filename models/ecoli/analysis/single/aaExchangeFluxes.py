from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Amino acid IDs
		sim_data = cPickle.load(open(simDataFile, "rb"))
		aaIDs = sim_data.moleculeGroups.aaIDs

		# Amino acid exchanges fluxes
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime
		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		externalExchangeFluxes = fba_results.readColumn("externalExchangeFluxes")
		externalMoleculeIDs = fba_results.readAttribute("externalMoleculeIDs")

		# Plot
		rows = 6
		cols = 4
		fig = plt.figure(figsize = (8, 11.5))

		for plotIndex, aa in enumerate(aaIDs):
			ax = plt.subplot(rows, cols, plotIndex + 1)

			if not aa.startswith("L-SELENOCYSTEINE"):
				aa = aa[:-3] + "[p]"
			if aa in externalMoleculeIDs:
				aaFlux = externalExchangeFluxes[:, externalMoleculeIDs.index(aa)]
			else:
				aaFlux = np.zeros(len(time))

			ax.plot(time / 60., aaFlux)
			ax.set_xlabel("Time (min)", fontsize = 6)
			ax.set_ylabel("mmol/gDCW/hr", fontsize = 6)
			ax.set_title("%s" % aa, fontsize = 6, y = 1.1)
			ax.tick_params(which = "both", direction = "out", labelsize = 6)

		plt.rc("font", size = 6)
		plt.suptitle("External exchange fluxes of amino acids", fontsize = 10)

		plt.subplots_adjust(hspace = 1, wspace = 1)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
