#!/usr/bin/env python
"""
Plots gene copy number

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/2/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "GeneCopyNumber.hdf"))

	table = h.root.GeneCopyNumber

	# geneCopyNumber = table.read(0, None, 1, 'gene_copy_number')
	totalCopyNumber = table.read(0, None, 1, 'total_copy_number')

	h.close()

	h = tables.open_file(os.path.join(simOutDir, "Mass.hdf"))
	table = h.root.Mass
	time = np.array([x["time"] for x in table.iterrows()])
	h.close()

	plt.figure(figsize = (8.5, 11))

	# plt.subplot(2,1,1)
	plt.plot(time / 60., totalCopyNumber, linewidth = 2)
	plt.xlabel("Time (min)")
	plt.ylabel("Count")
	plt.title("Total gene copy number")


	# plt.subplot(2,1,2)
	# plt.plot(time / 60., geneCopyNumber, linewidth = 2)
	# plt.xlabel("Time (min)")
	# plt.ylabel("Count")
	# plt.title("Individual gene copy number")

	# plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

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