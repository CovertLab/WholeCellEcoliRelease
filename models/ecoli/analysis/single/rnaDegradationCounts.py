#!/usr/bin/env python
"""
Plots counts of rna degraded and the resulting free NMPs

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/2015 - Updated 8/10/2015
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

import wholecell.utils.constants
from wholecell.io.tablereader import TableReader

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


def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))

	endoRnaseIds = kb.moleculeGroups.endoRnaseIds
	exoRnaseIds = kb.moleculeGroups.exoRnaseIds
	RnaseIds = np.concatenate((endoRnaseIds, exoRnaseIds))

	# Load count data for s30 proteins, rRNA, and final 30S complex
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")

	# Get indexes
	proteinIndexes = np.array([moleculeIds.index(protein) for protein in RnaseIds], np.int)
	exoproteinIndexes = np.array([moleculeIds.index(protein) for protein in exoRnaseIds], np.int)
	endoproteinIndexes = np.array([moleculeIds.index(protein) for protein in endoRnaseIds], np.int)

	# Load data
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	RnaseCounts = bulkMolecules.readColumn("counts")[:, proteinIndexes]

	exoRnaseCounts = bulkMolecules.readColumn("counts")[:, exoproteinIndexes]
	endoRnaseCounts = bulkMolecules.readColumn("counts")[:, endoproteinIndexes]
	bulkMolecules.close()

	rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
	countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
	nucleotidesFromDegradation = rnaDegradationListenerFile.readColumn('nucleotidesFromDegradation')
	FractionActiveEndoRNases = rnaDegradationListenerFile.readColumn('FractionActiveEndoRNases')
	rnaDegradationListenerFile.close()

	# Computation
	totalRnaseCounts = RnaseCounts.sum(axis = 1)
	requiredRnaseTurnover = nucleotidesFromDegradation / RnaseCounts.sum(axis = 1)

	totalexoRnaseCounts = exoRnaseCounts.sum(axis = 1)
	totalendoRnaseCounts = endoRnaseCounts.sum(axis = 1)
	kcatExoRnase = nucleotidesFromDegradation / exoRnaseCounts.sum(axis = 1)
	kcatEndoRnase = countRnaDegraded.sum(axis = 1) / endoRnaseCounts.sum(axis = 1)

	# Plotting
	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	countRnaDegraded_axis = plt.subplot(9,1,1)
	countRnaDegraded_axis.plot(time / 60., countRnaDegraded.sum(axis = 1))
	plt.ylabel("Counts of RNA degraded", fontsize = 7)
	plt.xlabel("Time (min)")

	nucleotidesFromDegradation_axis = plt.subplot(9,1,2)
	nucleotidesFromDegradation_axis.plot(time / 60., nucleotidesFromDegradation)
	plt.ylabel("Nucleotides from RNA degraded", fontsize = 7)
	plt.xlabel("Time (min)")

	totalRnaseCounts_axis = plt.subplot(9,1,3)
	totalRnaseCounts_axis.plot(time / 60., totalRnaseCounts)
	plt.ylabel("Total RNAse counts", fontsize = 7)
	plt.xlabel("Time (min)")

	requiredRnaseTurnover_axis = plt.subplot(9,1,4)
	requiredRnaseTurnover_axis.plot(time / 60., requiredRnaseTurnover)
	plt.ylabel("RNase turnover required (nt/s)", fontsize = 7)
	plt.xlabel("Time (min)")

	totalexoRnaseCounts_axis = plt.subplot(9,1,5)
	totalexoRnaseCounts_axis.plot(time / 60., totalexoRnaseCounts)
	plt.ylabel("Counts of exoRNase", fontsize = 7)
	plt.xlabel("Time (min)")

	totalendoRnaseCounts_axis = plt.subplot(9,1,6)
	totalendoRnaseCounts_axis.plot(time / 60., totalendoRnaseCounts)
	plt.ylabel("Counts of endoRNase", fontsize = 7)
	plt.xlabel("Time (min)")

	kcatExoRnase_axis = plt.subplot(9,1,7)
	kcatExoRnase_axis.plot(time / 60., kcatExoRnase)
	plt.ylabel("ExoRNase turnover required", fontsize = 7)
	plt.xlabel("Time (min)")

	kcatEndoRnase_axis = plt.subplot(9,1,8)
	kcatEndoRnase_axis.plot(time / 60., kcatEndoRnase)
	plt.ylabel("EndoRNase turnover required", fontsize = 7)
	plt.xlabel("Time (min)")

	FractionActiveEndoRNases_axis = plt.subplot(9,1,9)
	FractionActiveEndoRNases_axis.plot(time / 60., FractionActiveEndoRNases)
	plt.ylabel("Fraction of active EndoRNases", fontsize = 7)
	plt.xlabel("Time (min)")

	plt.subplots_adjust(hspace = 0.9, wspace = 0.5)

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
