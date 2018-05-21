#!/usr/bin/env python
"""
Plot rna synthesis probabilities

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/9/2016
"""

from __future__ import division

import argparse
import os

import numpy as np
import cPickle
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants


def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))
	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	isRRna = sim_data.process.transcription.rnaData["isRRna"]
	isTRna = sim_data.process.transcription.rnaData["isTRna"]

	rnaSynthProbListener = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
	rnaIds = rnaSynthProbListener.readAttribute('rnaIds')
	rnaSynthProb = rnaSynthProbListener.readColumn('rnaSynthProb')
	time = rnaSynthProbListener.readColumn('time')
	rnaSynthProbListener.close()

	mRnaSynthProb = rnaSynthProb[:, isMRna].sum(axis = 1)
	rRnaSynthProb = rnaSynthProb[:, isRRna].sum(axis = 1)
	tRnaSynthProb = rnaSynthProb[:, isTRna].sum(axis = 1)



	# Plot
	rows = 3
	cols = 1
	fig = plt.figure(figsize = (11, 8.5))
	plt.figtext(0.4, 0.96, "RNA synthesis probabilities over time", fontsize = 12)
	nMRnas = np.sum(isMRna)
	nRRnas = np.sum(isRRna)
	nTRnas = np.sum(isTRna)
	subplotOrder = [mRnaSynthProb, rRnaSynthProb, tRnaSynthProb]
	subplotTitles = ["mRNA\n(sum of %s mRNAs)" % nMRnas, "rRNA\n(sum of %s rRNAs)" % nRRnas, "tRNA\n(sum of %s tRNAs)" % nTRnas]

	for index, rnaSynthProb in enumerate(subplotOrder):
		ax = plt.subplot(rows, cols, index + 1)
		ax.plot(time, rnaSynthProb)

		ax.set_title(subplotTitles[index], fontsize = 10)
		ymin = np.min(rnaSynthProb)
		ymax = np.max(rnaSynthProb)
		yaxisBuffer = np.around(1.2*(ymax - ymin), 3)
		ax.set_ylim([ymin, yaxisBuffer])
		ax.set_yticks([ymin, ymax, yaxisBuffer])
		ax.set_yticklabels([ymin, np.around(ymax, 3), yaxisBuffer], fontsize = 10)
		ax.set_xlim([time[0], time[-1]])
		ax.tick_params(axis = "x", labelsize = 10)
		ax.spines["left"].set_visible(False)
		ax.spines["right"].set_visible(False)


	plt.subplots_adjust(hspace = 0.5, )
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	plt.close("all")


if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)
	parser.add_argument("--validationDataFile", help = "", type = str, )

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
