from __future__ import absolute_import

import os

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

		# Exchange flux
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		exFlux = fba_results.readColumn("externalExchangeFluxes")
		exMolec = fba_results.readAttribute("externalMoleculeIDs")
		moleculeIDs = ["GLC[p]", "OXYGEN-MOLECULE[p]"]

		# Plot
		fig = plt.figure(figsize = (8, 11.5))
		rows = len(moleculeIDs)
		cols = 1

		for index, molecule in enumerate(["GLC[p]", "OXYGEN-MOLECULE[p]"]):
			if molecule not in exMolec:
				continue
			moleculeFlux = -1. * exFlux[:, exMolec.index(molecule)]
			ax = plt.subplot(rows, cols, index + 1)
			ax.plot(time / 60. / 60., moleculeFlux)

			averageFlux = np.average(moleculeFlux)
			yRange = np.min([np.abs(np.max(moleculeFlux) - averageFlux), np.abs(np.min(moleculeFlux) - averageFlux)])
			ymin = np.round(averageFlux - yRange)
			ymax = np.round(averageFlux + yRange)
			ax.set_ylim([ymin, ymax])

			abs_max = np.max(moleculeFlux)
			abs_min = np.min(moleculeFlux)

			plt.figtext(0.7, 1. / float(rows) * 0.7 + (rows - 1 - index) / float(rows),
				"Max: %s\nMin: %s" % (abs_max, abs_min), fontsize = 8)

			ax.set_ylabel("External %s\n(mmol/gDCW/hr)" % molecule, fontsize = 8)
			ax.set_xlabel("Time (hr)", fontsize = 8)
			ax.set_title("%s" % molecule, fontsize = 10, y = 1.1)
			ax.tick_params(labelsize = 8, which = "both", direction = "out")


		plt.subplots_adjust(hspace = 0.5, wspace = 1)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
