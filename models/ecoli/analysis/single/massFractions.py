"""
Plot mass fractions

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import absolute_import

import os

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

		mass = TableReader(os.path.join(simOutDir, "Mass"))

		cell = mass.readColumn("cellMass")
		cellDry = mass.readColumn("dryMass")
		protein = mass.readColumn("proteinMass")
		rna = mass.readColumn("rnaMass")
		smallMolecules = mass.readColumn("smallMoleculeMass")
		# tRna = mass.readColumn("tRnaMass")
		# rRna = mass.readColumn("rRnaMass")
		# mRna = mass.readColumn("mRnaMass")
		# dna = mass.readColumn("dnaMass")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		t = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime


		mass.close()

		plt.figure(figsize = (8.5, 11))

		plt.subplot(5, 1, 1)

		plt.plot(t / 60., cell, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * cell[0], 2 * cell[0]], 'r--')
		plt.xlabel("Time (min)")
		plt.ylabel("Total Mass (fg)")
		plt.title("Total Mass Final:Initial = %0.2f" % (cell[-1] / cell[0]))

		plt.subplot(5, 1, 2)

		plt.plot(t / 60., cellDry, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * cellDry[0], 2 * cellDry[0]], 'r--')
		plt.xlabel("Time (min)")
		plt.ylabel("Dry Mass (fg)")
		plt.title("Dry Mass Final:Initial = %0.2f" % (cellDry[-1] / cellDry[0]))

		plt.subplot(5, 1, 3)

		plt.plot(t / 60., protein, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * protein[0], 2 * protein[0]], "r--")
		plt.xlabel("Time (min)")
		plt.ylabel("Protein Mass (fg)")
		plt.title("Total Protein Mass Final:Initial = %0.2f" % (protein[-1] / protein[0]))
		plt.show()

		plt.subplot(5, 1, 4)

		plt.plot(t / 60., rna, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * rna[0], 2 * rna[0]], "r--")
		plt.xlabel("Time (min)")
		plt.ylabel("RNA Mass (fg)")
		plt.title("Total RNA Mass Final:Initial = %0.2f" % (rna[-1] / rna[0]))

		plt.subplot(5, 1, 5)

		plt.plot(t / 60., smallMolecules, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * smallMolecules[0], 2 * smallMolecules[0]], "r--")
		plt.xlabel("Time (min)")
		plt.ylabel("Small molecules (fg)")
		plt.title("Total Small Molecule Mass Final:Initial = %0.2f" % (smallMolecules[-1] / smallMolecules[0]))

		plt.subplots_adjust(hspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
