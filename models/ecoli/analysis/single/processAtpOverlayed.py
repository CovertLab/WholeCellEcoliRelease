#!/usr/bin/env python
"""
Plot water allocation for each process

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/20/2016
"""

from __future__ import division

import argparse
import os

import numpy as np
from scipy.stats import pearsonr
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils.fitting import normalize
from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.analysis.analysis_tools import exportFigure

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	processNames = bulkMolecules.readAttribute("processNames")

	atpAllocatedInitial = bulkMolecules.readColumn("atpAllocatedInitial")
	atpRequested = bulkMolecules.readColumn("atpRequested")

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	bulkMolecules.close()


	# Plot
	plt.figure(figsize = (8.5, 11))
	rows = 7
	cols = 2

	for processIndex in np.arange(len(processNames)):
		ax = plt.subplot(rows, cols, processIndex + 1)
		ax.plot(time / 60., atpAllocatedInitial[:, processIndex])
		ax.plot(time / 60., atpRequested[:, processIndex])
		ax.set_title(str(processNames[processIndex]), fontsize = 8, y = 0.85)

		ymin = np.amin([atpAllocatedInitial[:, processIndex], atpRequested[:, processIndex]])
		ymax = np.amax([atpAllocatedInitial[:, processIndex], atpRequested[:, processIndex]])
		ax.set_ylim([ymin, ymax])
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		# ax.set_xticks([])

	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

	plt.subplots_adjust(hspace = 2.0, wspace = 2.0)

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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
