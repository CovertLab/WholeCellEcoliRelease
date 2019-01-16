from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
	[228,26,28],
	[55,126,184],
	[77,175,74],
	[152,78,163],
	[255,127,0],
	[255,255,51],
	[166,86,40],
	[247,129,191]
	]

COLORS = [
	[colorValue/255. for colorValue in color]
	for color in COLORS_256
	]


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		mass = TableReader(os.path.join(simOutDir, "Mass"))

		# cell = mass.readColumn("cellMass")
		# cellDry = mass.readColumn("dryMass")
		protein = mass.readColumn("proteinMass")
		# rna = mass.readColumn("rnaMass")
		tRna = mass.readColumn("tRnaMass")
		rRna = mass.readColumn("rRnaMass")
		mRna = mass.readColumn("mRnaMass")
		dna = mass.readColumn("dnaMass")
		smallMolecules = mass.readColumn("smallMoleculeMass")

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		t = main_reader.readColumn("time") - initialTime


		masses = np.vstack([
			protein/protein[0],
			rRna/rRna[0],
			tRna/tRna[0],
			mRna/mRna[0],
			dna/dna[0],
			smallMolecules/smallMolecules[0],
			]).T

		massLabels = ["Protein", "rRNA", "tRNA", "mRNA", "DNA", "Small Mol.s"]

		plt.figure(figsize = (8.5, 11))

		# plt.rc('axes', color_cycle=COLORS)
		plt.gca().set_prop_cycle('color', COLORS)

		plt.plot(t / 60., masses, linewidth = 2)
		plt.xlabel("Time (min)")
		plt.ylabel("Mass (normalized by t = 0 min)")
		plt.title("Biomass components")
		#plt.axis([0, 60, 0.5, 2.5])

		plt.legend(massLabels, loc = "best")

		# plt.show()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
