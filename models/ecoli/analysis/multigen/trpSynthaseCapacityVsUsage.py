#!/usr/bin/env python
"""
Plot trp regulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/2016
"""

from __future__ import division

import argparse
import os

import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils.fitting import normalize
from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths

BURN_IN = 10

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	allDirs = ap.get_cells()

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))
	trpIdx = sim_data.moleculeGroups.aaIDs.index("TRP[c]")

	plt.figure(figsize = (8.5, 11))

	for simDir in allDirs:
		simOutDir = os.path.join(simDir, "simOut")

		growthLimits = TableReader(os.path.join(simOutDir, "GrowthLimits"))

		trpRequests = growthLimits.readColumn("aaRequestSize")[BURN_IN:, trpIdx]

		growthLimits.close()

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

		moleculeIds = bulkMolecules.readAttribute("objectNames")

		trpSynIdx = moleculeIds.index("TRYPSYN[c]")

		trpSynCounts = bulkMolecules.readColumn("counts")[BURN_IN:, trpSynIdx]

		bulkMolecules.close()

		trpSynKcat = 2**( (37. - 25.) / 10.) * 4.1 # From PMID 6402362 (kcat of 4.1/s measured at 25 C)

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")[BURN_IN:]
		timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")[BURN_IN:]


		trpSynMaxCapacity = trpSynKcat * trpSynCounts * timeStep

		
		##############################################################
		ax = plt.subplot(3, 1, 1)
		ax.plot(time / 60., trpSynMaxCapacity, color = "b")
		plt.ylabel("Tryptophan Synthase Max Capacity", fontsize = 10)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 10)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(3, 1, 2)
		ax.plot(time, trpRequests, color = "b")
		plt.ylabel("Trp Requested By Translation", fontsize = 10)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize=10)
		ax.set_xticks([])
		##############################################################


		##############################################################
		ax = plt.subplot(3, 1, 3)
		ax.plot(time / 3600., trpSynMaxCapacity / trpRequests, color = "b")
		ax.plot([0, time[-1] / 3600.], [1., 1.], "k--")
		plt.ylabel("(Max capacity) / (Request)", fontsize = 10)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.2f" % ymin, "%0.2f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		# ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 10)
		ax.set_xticks(ax.get_xlim())
		##############################################################

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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
