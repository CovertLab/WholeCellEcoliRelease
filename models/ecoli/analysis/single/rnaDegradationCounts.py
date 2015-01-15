#!/usr/bin/env python
"""
Plots counts of rna degraded and the resulting free NMPs

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/2015
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

import wholecell.utils.constants

FONT = {
		'size'	:	8
		}

def setAxisMaxMin(axis, data):
	ymax = np.max(data)
	ymin = 0
	if ymin == ymax:
		axis.set_yticks([ymin])
	else:
		axis.set_yticks([ymin, ymax])

def sparklineAxis(axis, x, y, tickPos, lineType, color):
	axis.plot(x, y, linestyle = 'steps' + lineType, color = color, linewidth = 2)
	axis.spines['top'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.yaxis.set_ticks_position(tickPos)
	axis.xaxis.set_ticks_position('none')
	axis.tick_params(which = 'both', direction = 'out')
	axis.tick_params(labelbottom = 'off')
	for tl in axis.get_yticklabels():
		tl.set_color(color)


def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	RnaseIds = ["EG10856-MONOMER[p]", "EG11620-MONOMER[c]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]",
	"EG10858-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]",
	"G7365-MONOMER[c]", "EG10862-MONOMER[c]", "EG10863-MONOMER[c]", "EG11259-MONOMER[c]", "EG11547-MONOMER[c]"]
	

	# Load count data for s30 proteins, rRNA, and final 30S complex
	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:
		# Get indexes
		moleculeIds = bulkMoleculesFile.root.names.moleculeIDs.read()
		proteinIndexes = np.array([moleculeIds.index(protein) for protein in RnaseIds], np.int)

		# Load data
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules
		time = bulkMolecules.col("time")
		RnaseCounts = bulkMolecules.read(0, None, 1, "counts")[:, proteinIndexes]

	with tables.open_file(os.path.join(simOutDir, "RnaDegradationListener.hdf")) as rnaDegradationListenerFile:
		time = rnaDegradationListenerFile.root.RnaDegradationListener.col('time')
		countRnaDegraded = rnaDegradationListenerFile.root.RnaDegradationListener.col('countRnaDegraded')
		nucleotidesFromDegradation = rnaDegradationListenerFile.root.RnaDegradationListener.col('nucleotidesFromDegradation')

	# Computation
	totalRnaseCounts = RnaseCounts.sum(axis = 1)
	requiredRnaseTurnover = nucleotidesFromDegradation / RnaseCounts.sum(axis = 1)

	# Plotting
	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	countRnaDegraded_axis = plt.subplot(4,1,1)
	countRnaDegraded_axis.plot(time / 60, countRnaDegraded)

	nucleotidesFromDegradation_axis = plt.subplot(4,1,2)
	nucleotidesFromDegradation_axis.plot(time / 60, nucleotidesFromDegradation)

	totalRnaseCounts_axis = plt.subplot(4,1,3)
	totalRnaseCounts_axis.plot(time / 60, totalRnaseCounts)

	requiredRnaseTurnover_axis = plt.subplot(4,1,4)
	requiredRnaseTurnover_axis.plot(time / 60, requiredRnaseTurnover)

	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

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
