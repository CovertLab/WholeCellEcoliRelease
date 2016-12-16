#!/usr/bin/env python
"""
Plot ATP allocation for each process

@author: Heejo Choi
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
from wholecell.analysis.analysis_tools import exportFigure

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))
	nAvogadro = sim_data.constants.nAvogadro
	cellDensity = sim_data.constants.cellDensity
	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = units.fg * mass.readColumn("cellMass")
	cellVolume = cellMass / cellDensity

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	bulkMoleculeIds = bulkMolecules.readAttribute("objectNames")

	processNames = bulkMolecules.readAttribute("processNames")
	
	atpAllocatedInitial = bulkMolecules.readColumn("atpAllocatedInitial")
	atpRequested = bulkMolecules.readColumn("atpRequested")

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	atpId = np.where([str(x) == "ATP[c]" for x in bulkMoleculeIds])[0][0]
	atpCounts = bulkMolecules.readColumn("counts")[:, atpId]
	atpConcentrations = atpCounts / (cellVolume.asNumber(units.L) * nAvogadro.asNumber(1 / units.mol))

	bulkMolecules.close()

	# Plot
	plt.figure(figsize = (8.5, 11))
	rows = 14
	cols = 2
	processIndex = 0

	for subplotIndex in np.arange(0, 2* len(processNames) -1, 2):
		ax = plt.subplot(rows, cols, subplotIndex + 1)
		ax.plot(time / 60., atpAllocatedInitial[:, processIndex], color = "b")
		ax.set_title(str(processNames[processIndex]) + " Allocated", fontsize = 8, y = 0.80)
		ymin = np.amin(atpAllocatedInitial[:, processIndex])
		ymax = np.amax(atpAllocatedInitial[:, processIndex])
		ax.set_ylim([ymin, ymax])
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])

		ax.axvline(x = time[916] / 60., color = "r", linestyle = "--")


		ax = plt.subplot(rows, cols, subplotIndex + 2)
		ax.plot(time / 60., atpRequested[:, processIndex], color = "g")
		ax.set_title(str(processNames[processIndex]) + " Requested", fontsize = 8, y = 0.80)
		ymin = np.amin(atpRequested[:, processIndex])
		ymax = np.amax(atpRequested[:, processIndex])
		ax.set_ylim([ymin, ymax])
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		ax.axvline(x = time[916] / 60., color = "r", linestyle = "--")

		processIndex += 1

	# Plot ATP counts
	ax = plt.subplot(rows, cols, 27)
	ax.plot(time / 60., atpCounts)
	ax.set_title("ATP[c] counts", fontsize = 8)
	ymin = np.amin(atpCounts)
	ymax = np.amax(atpCounts)
	ax.set_ylim([ymin, ymax])
	ax.set_yticks([ymin, ymax])
	ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
	# ax.set_xticks([])
	ax.axvline(x = time[916] / 60., color = "r", linestyle = "--")

	# Plot ATP concentrations
	ax = plt.subplot(rows, cols, 28)
	ax.plot(time / 60., atpConcentrations)
	ax.set_title("ATP[c] concentrations", fontsize = 8)
	ymin = np.amin(atpConcentrations)
	ymax = np.amax(atpConcentrations)
	ax.set_ylim([ymin, ymax])
	ax.set_yticks([ymin, ymax])
	ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
	# ax.set_xticks([])
	ax.axvline(x = time[916] / 60., color = "r", linestyle = "--")


	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

	plt.subplots_adjust(hspace = 1, wspace = 1, bottom = 0.2, top = 0.8, left = 0.2, right = 0.8)

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
