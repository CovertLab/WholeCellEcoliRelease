#!/usr/bin/env python
"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/27/2014
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

THRESHOLD = 1e-13 # roughly, the mass of an electron

FG_PER_DALTON = 1.6605402e-9

# TODO: get these from the KB
REPRESENTATIVE_MASSES = {
	"proton":1.007 * FG_PER_DALTON,
	"amino acid":109 * FG_PER_DALTON,
	"ATP":551 * FG_PER_DALTON,
	"protein":40e3 * FG_PER_DALTON,
	"ribosome":2700e3 * FG_PER_DALTON
	}

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	mass = TableReader(os.path.join(simOutDir, "Mass"))

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	processMassDifferences = mass.readColumn("processMassDifferences")

	processNames = mass.readAttribute("processNames")

	mass.close()

	avgProcessMassDifferences = np.abs(processMassDifferences).sum(axis = 0) / len(time)

	index = np.arange(len(processNames))
	width = 1

	plt.figure(figsize = (8.5, 11))

	axes = plt.axes()

	r1 = axes.barh(index, avgProcessMassDifferences * (avgProcessMassDifferences > THRESHOLD), width, log = True, color = (0.9, 0.2, 0.2))
	r2 = axes.barh(index, avgProcessMassDifferences * (avgProcessMassDifferences <= THRESHOLD), width, log = True, color = (0.2, 0.2, 0.9))

	axes.set_yticks(index+width/2)
	axes.set_yticklabels(processNames) #, rotation = -45)

	axes.plot([THRESHOLD, THRESHOLD], [index[0], index[-1]+width], 'k--', linewidth=3)

	plt.text(THRESHOLD, index[-1], "electron", rotation = "vertical", va = "center", ha = "right")

	for name, mass in REPRESENTATIVE_MASSES.viewitems():
		plt.axvline(mass, color = "k")
		plt.text(mass, index[-1], name, rotation = "vertical", va = "center", ha = "right")

	plt.xlabel("Mass difference (fg)")

	plt.title("Average absolute change in mass by individual processes")

	plt.tight_layout()
	plt.grid(True, which = "major")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")


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
