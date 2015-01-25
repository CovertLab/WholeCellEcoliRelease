#!/usr/bin/env python
"""
Plot mass fractions

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	mass = TableReader(os.path.join(simOutDir, "Mass"))

	cell = mass.readColumn("cellMass")
	cellDry = mass.readColumn("dryMass")
	protein = mass.readColumn("proteinMass")
	rna = mass.readColumn("rnaMass")
	# tRna = mass.readColumn("tRnaMass")
	# rRna = mass.readColumn("rRnaMass")
	# mRna = mass.readColumn("mRnaMass")
	# dna = mass.readColumn("dnaMass")
	t = mass.readColumn("time")

	mass.close()

	plt.figure(figsize = (8.5, 11))

	plt.subplot(4, 1, 1)

	plt.plot(t / 60., cell, linewidth = 2)
	plt.plot([t[0] / 60., t[-1] / 60.], [2 * cell[0], 2 * cell[0]], 'r--')
	plt.xlabel("Time (min)")
	plt.ylabel("Total Mass (fg)")
	plt.title("Total Mass Final:Initial = %0.2f" % (cell[-1] / cell[0]))

	plt.subplot(4, 1, 2)

	plt.plot(t / 60., cellDry, linewidth = 2)
	plt.plot([t[0] / 60., t[-1] / 60.], [2 * cellDry[0], 2 * cellDry[0]], 'r--')
	plt.xlabel("Time (min)")
	plt.ylabel("Dry Mass (fg)")
	plt.title("Dry Mass Final:Initial = %0.2f" % (cellDry[-1] / cellDry[0]))

	plt.subplot(4, 1, 3)

	plt.plot(t / 60., protein, linewidth = 2)
	plt.plot([t[0] / 60., t[-1] / 60.], [2 * protein[0], 2 * protein[0]], "r--")
	plt.xlabel("Time (min)")
	plt.ylabel("Protein Mass (fg)")
	plt.title("Total Protein Mass Final:Initial = %0.2f" % (protein[-1] / protein[0]))
	plt.show()

	plt.subplot(4, 1, 4)

	plt.plot(t / 60., rna, linewidth = 2)
	plt.plot([t[0] / 60., t[-1] / 60.], [2 * rna[0], 2 * rna[0]], "r--")
	plt.xlabel("Time (min)")
	plt.ylabel("RNA Mass (fg)")
	plt.title("Total RNA Mass Final:Initial = %0.2f" % (rna[-1] / rna[0]))

	plt.subplots_adjust(hspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)


if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
