"""
Plot mass fractions

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		mass = TableReader(os.path.join(simOutDir, "Mass"))

		cell = mass.readColumn("cellMass")
		cellDry = mass.readColumn("dryMass")
		protein = mass.readColumn("proteinMass")
		rna = mass.readColumn("rnaMass")
		smallMolecules = mass.readColumn("smallMoleculeMass")
		dna = mass.readColumn("dnaMass")

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		t = main_reader.readColumn("time") - initialTime

		fig = plt.figure(figsize = (8.5, 15))
		n_subplots = 6
		second_axis_color = 'g'

		plt.subplot(n_subplots, 1, 1)
		plt.plot(t / 60., cell, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * cell[0], 2 * cell[0]], 'r--')
		plt.ylabel("Total Mass (fg)")
		plt.title("Total Mass Final:Initial = %0.2f" % (cell[-1] / cell[0]), fontsize=8)

		plt.subplot(n_subplots, 1, 2)
		plt.plot(t / 60., cellDry, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * cellDry[0], 2 * cellDry[0]], 'r--')
		plt.ylabel("Dry Mass (fg)")
		plt.title("Dry Mass Final:Initial = %0.2f" % (cellDry[-1] / cellDry[0]), fontsize=8)

		ax = plt.subplot(n_subplots, 1, 3)
		plt.plot(t / 60., protein, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * protein[0], 2 * protein[0]], "r--")
		plt.ylabel("Protein Mass (fg)")
		plt.title("Total Protein Mass Final:Initial = %0.2f\nAverage dry mass fraction: %0.3f"
			% (protein[-1] / protein[0], np.mean(protein / cellDry)), fontsize=8)
		ax2 = ax.twinx()
		ax2.plot(t / 60., protein / cellDry, second_axis_color)
		ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
		ax2.set_yticks(ax2.get_ylim())
		ax2.set_ylabel('Fraction of dry mass', color=second_axis_color)

		ax = plt.subplot(n_subplots, 1, 4)
		plt.plot(t / 60., rna, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * rna[0], 2 * rna[0]], "r--")
		plt.ylabel("RNA Mass (fg)")
		plt.title("Total RNA Mass Final:Initial = %0.2f\nAverage dry mass fraction: %0.3f"
			% (rna[-1] / rna[0], np.mean(rna / cellDry)), fontsize=8)
		ax2 = ax.twinx()
		ax2.plot(t / 60., rna / cellDry, second_axis_color)
		ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
		ax2.set_yticks(ax2.get_ylim())
		ax2.set_ylabel('Fraction of dry mass', color=second_axis_color)

		ax = plt.subplot(n_subplots, 1, 5)
		plt.plot(t / 60., smallMolecules, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * smallMolecules[0], 2 * smallMolecules[0]], "r--")
		plt.ylabel("Small molecules (fg)")
		plt.title("Total Small Molecule Mass Final:Initial = %0.2f\nAverage dry mass fraction: %0.3f"
			% (smallMolecules[-1] / smallMolecules[0], np.mean(smallMolecules / cellDry)), fontsize=8)
		ax2 = ax.twinx()
		ax2.plot(t / 60., smallMolecules / cellDry, second_axis_color)
		ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
		ax2.set_yticks(ax2.get_ylim())
		ax2.set_ylabel('Fraction of dry mass', color=second_axis_color)

		ax = plt.subplot(n_subplots, 1, 6)
		plt.plot(t / 60., dna, linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [2 * dna[0], 2 * dna[0]], "r--")
		plt.xlabel("Time (min)")
		plt.ylabel("DNA Mass (fg)")
		plt.title("Total DNA Mass Final:Initial = %0.2f\nAverage dry mass fraction: %0.3f"
			% (dna[-1] / dna[0], np.mean(dna / cellDry)), fontsize=8)
		ax2 = ax.twinx()
		ax2.plot(t / 60., dna / cellDry, second_axis_color)
		ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
		ax2.set_yticks(ax2.get_ylim())
		ax2.set_ylabel('Fraction of dry mass', color=second_axis_color)

		plt.subplots_adjust(hspace = 0.75)
		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
