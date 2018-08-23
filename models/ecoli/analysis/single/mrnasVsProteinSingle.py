from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

# TODO: account for complexation

CMAP_COLORS_255 = [
	[103,0,31],
	[178,24,43],
	[214,96,77],
	[244,165,130],
	[253,219,199],
	[247,247,247],
	[209,229,240],
	[146,197,222],
	[67,147,195],
	[33,102,172],
	[5,48,97],
	]

CMAP_COLORS = [[shade/255. for shade in color] for color in CMAP_COLORS_255]


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get the names of rnas from the KB

		sim_data = cPickle.load(open(simDataFile, "rb"))

		rnaIds = sim_data.process.transcription.rnaData["id"][sim_data.relation.rnaIndexToMonomerMapping]

		proteinIds = sim_data.process.translation.monomerData["id"]

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")

		moleculeIds = bulkMolecules.readAttribute("objectNames")

		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

		rnaCountsBulk = bulkMoleculeCounts[:, rnaIndexes]

		proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)

		proteinCountsBulk = bulkMoleculeCounts[:, proteinIndexes]

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		bulkMolecules.close()

		time = time/60.

		# TODO: plots for all molecules (needs to be fast)
		# TODO: account for complex composition, other modifications

		index = 611
		name = "DNA polymerase"

		plt.figure(figsize = (8.5, 11))

		# TODO: functionalize this as a filled step plot

		stepTime = np.roll(np.repeat(time, 2), +1)[1:]

		stepRna = np.repeat(rnaCountsBulk[:, index], 2)[1:]
		stepProtein = np.repeat(proteinCountsBulk[:, index], 2)[1:]

		plt.subplot(2, 1, 1)

		polygonRna = np.vstack([
			np.hstack([stepTime[0], stepTime, stepTime[-1]]),
			np.hstack([0, stepRna, 0])
			]).T

		plt.gca().add_patch(plt.Polygon(
			polygonRna,
			facecolor = CMAP_COLORS[3],
			edgecolor = "none"
			))

		lineRna = np.vstack([stepTime, stepRna]).T

		plt.gca().add_patch(plt.Polygon(
			lineRna,
			facecolor = "none",
			edgecolor = CMAP_COLORS[1],
			closed = False,
			linewidth = 2
			))

		plt.axis("auto")

		plt.xlabel("Time (min)")
		plt.ylabel("mRNA copies (number)")

		plt.subplot(2, 1, 2)

		polygonProtein = np.vstack([
			np.hstack([stepTime[0], stepTime, stepTime[-1]]),
			np.hstack([0, stepProtein, 0])
			]).T

		plt.gca().add_patch(plt.Polygon(
			polygonProtein,
			facecolor = CMAP_COLORS[-4],
			edgecolor = "none"
			))

		lineProtein = np.vstack([stepTime, stepProtein]).T

		plt.gca().add_patch(plt.Polygon(
			lineProtein,
			facecolor = "none",
			edgecolor = CMAP_COLORS[-2],
			closed = False,
			linewidth = 2
			))

		plt.axis("auto")

		plt.xlabel("Time (min)")
		plt.ylabel("Protein copies (number)")

		plt.suptitle(name)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
