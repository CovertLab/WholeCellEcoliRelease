from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils.sparkline import whitePadSparklineAxis


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		mass = TableReader(os.path.join(simOutDir, "Mass"))

		protein = mass.readColumn("proteinMass")
		tRna = mass.readColumn("tRnaMass")
		rRna = mass.readColumn("rRnaMass")
		mRna = mass.readColumn("mRnaMass")
		dna = mass.readColumn("dnaMass")
		smallMolecules = mass.readColumn("smallMoleculeMass")

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		t = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		masses = np.vstack([
			protein/protein[0],
			rRna/rRna[0],
			tRna/tRna[0],
			mRna/mRna[0],
			dna/dna[0],
			smallMolecules/smallMolecules[0],
			]).T

		massLabels = ["Protein", "rRNA", "tRNA", "mRNA", "DNA", "Small Mol."]

		plt.figure(figsize = (2, 2.5))

		ax = plt.gca()
		ax.set_prop_cycle(plt.style.library['fivethirtyeight']['axes.prop_cycle'])

		plt.plot(t / 60., masses, linewidth=1)
		plt.xlabel("Time (min)")
		plt.ylabel("Mass (normalized by t = 0 min)")
		plt.title("Biomass components")
		plt.legend(massLabels, loc = "best")
		plt.axhline(2, linestyle='--', color='k')

		whitePadSparklineAxis(ax)

		xticks = [0, t[-1] / 60]
		yticks = [1, 2, 2.2]
		ax.set_xlim(xticks)
		ax.set_ylim((yticks[0], yticks[-1]))
		ax.set_xticks(xticks)
		ax.set_yticks(yticks)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		ax.set_xlabel("")
		ax.set_ylabel("")
		ax.set_title("")
		ax.set_xticklabels([])
		ax.set_yticklabels([])
		ax.legend_.remove()

		exportFigure(plt, plotOutDir, plotOutFileName + "_stripped", metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
