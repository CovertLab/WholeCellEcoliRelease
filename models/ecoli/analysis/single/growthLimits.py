#!/usr/bin/env python
"""
Plots various effects that may be limiting growth

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/2015
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	nAvogadro = kb.constants.nAvogadro
	moleculeIds = kb.moleculeGroups.aaIDs
	moleculeIds.append('GTP[c]')

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Load data
	growthLimitsDataFile = TableReader(os.path.join(simOutDir, "GrowthLimits"))

	gtpPoolSize = growthLimitsDataFile.readColumn("gtpPoolSize")
	gtpRequestSize = growthLimitsDataFile.readColumn("gtpRequestSize")
	gtpAllocated = growthLimitsDataFile.readColumn("gtpAllocated")
	gtpUsed = growthLimitsDataFile.readColumn("gtpUsed")

	aaPoolSize = growthLimitsDataFile.readColumn("aaPoolSize")
	aaRequestSize = growthLimitsDataFile.readColumn("aaRequestSize")
	aaAllocated = growthLimitsDataFile.readColumn("aaAllocated")
	aasUsed = growthLimitsDataFile.readColumn("aasUsed")

	poolSize = np.hstack((aaPoolSize,gtpPoolSize.reshape(gtpPoolSize.size,1)))
	requestSize = np.hstack((aaRequestSize,gtpRequestSize.reshape(gtpPoolSize.size,1)))
	allocated = np.hstack((aaAllocated,gtpAllocated.reshape(gtpPoolSize.size,1)))
	used =  np.hstack((aasUsed,gtpUsed.reshape(gtpPoolSize.size,1)))

	growthLimitsDataFile.close()

	# bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	# bulkMoleculeIds = bulkMolecules.readAttribute("objectNames")
	# bulkMolecules.close()

	# Plot
	fig = plt.figure(figsize = (11, 11))

	for idx in range(len(moleculeIds)):
		ax = plt.subplot(6,4,idx+1)
		ax.plot(time / 60., poolSize[:,idx], linewidth=2, label="pool size", color='k')
		ax.plot(time / 60., requestSize[:,idx], linewidth=2, label="request size", color='b')
		ax.plot(time / 60., allocated[:,idx], linewidth=2, label="allocated", color='r', linestyle='--')
		ax.plot(time / 60., used[:, idx], linewidth=2, label="used", color="c", linestyle=":")

		# Highlight title if request is greater than pool
		bbox = None
		if (poolSize < requestSize)[:,idx].any() or (allocated < requestSize)[:,idx].any():
			bbox = {'facecolor':'red', 'alpha':0.5, 'pad':10}
		elif (used[:,idx] < allocated[:,idx]).any():
			bbox = {'facecolor':'orange', 'alpha':0.5, 'pad':10}
		ax.set_title(moleculeIds[idx], bbox=bbox)

		# Set ticks so this is all easy to read
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize=10)
		ax.set_xticks([time.min() / 60., time.max() / 60.])
		ax.set_yticks([0., np.max([poolSize[:,idx].max(), requestSize[:,idx].max(), allocated[:,idx].max()])])

	# Create legend
	ax = plt.subplot(6,4,len(moleculeIds) + 1)
	ax.plot(0, 0, linewidth=2, label="pool size", color='k')
	ax.plot(0, 0, linewidth=2, label="request size", color='b')
	ax.plot(0, 0, linewidth=2, label="allocated", color='r', linestyle='--')
	ax.plot(0, 0, linewidth=2, label="used", color="c", linestyle=':')
	ax.legend(loc = 10,prop={'size':10})
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_title("Highlights process under-usage", fontsize=12, bbox={'facecolor':'orange', 'alpha':0.5, 'pad':10})
	ax.set_xlabel("Highlights pool/allocation limit", fontsize=12, bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

	# ax = plt.subplot(6,1,1)
	# ax.plot(time / 60., fractionGtpLimit, linewidth=2, color='b')
	# ax.set_ylabel("gtp used /\n gtp allocated")
	# ax.set_title("Fraction of GTP limit used by translation")

	# ax = plt.subplot(6,1,2)
	# ax.plot(time / 60., fractionAAsUsed, linewidth=2, color='b')
	# ax.set_ylabel("aa used /\n aa allocated")
	# ax.set_title("Fraction of AA allocated used by translation")

	# ax = plt.subplot(6,1,3)
	# ax.plot(time / 60., gtpPoolSize, linewidth=2, label="pool size", color='k')
	# ax.plot(time / 60., gtpRequestSize, linewidth=2, label="request size", color='b')
	# ax.plot(time / 60., gtpAllocated, linewidth=2, label="allocated", color='r', linestyle='--')
	# ax.legend()
	# ax.set_ylabel("GTP")

	# Save
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

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
