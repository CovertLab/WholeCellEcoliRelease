#!/usr/bin/env python
"""
Plot RNA polymerase terminations for fast transcription (rRNA and tRNA)

@date: Created 6/24/2015
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
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

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	rnapData = TableReader(os.path.join(simOutDir, "RnapData"))

	nTranscriptsTerminatedFast = rnapData.readColumn("nTerminatedFast")
	nTranscriptsTerminatedSlow = rnapData.readColumn("nTerminatedSlow")

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	rnapData.close()

	plt.figure(figsize = (8.5, 11))
	plt.title("mRNA versus t/rRNA terminations")

	plt.subplot(3,1,1)

	plt.plot(time / 60, nTranscriptsTerminatedFast,'r')
	plt.plot(time / 60, nTranscriptsTerminatedSlow, 'b')

	plt.xlabel("Time (min)")
	plt.ylabel("Number of RNA trascripts terminating")
	plt.legend("t or rRNA", "mRNA")

	plt.subplot(3,1,2)
	plt.plot(time / 60, nTranscriptsTerminatedFast, 'r')

	plt.xlabel("Time (min)")
	plt.ylabel("Number of t/rRNA trascripts terminating")

	plt.subplot(3,1,3)
	plt.plot(time / 60, nTranscriptsTerminatedSlow)

	plt.xlabel("Time (min)")
	plt.ylabel("Number of mRNA trascripts terminating")


	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

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

	import ipdb; ipdb.set_trace()
	
	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
