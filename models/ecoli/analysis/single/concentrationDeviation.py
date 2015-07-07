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

COLOR_CHOICES = np.array([
[0,0,0],
[252,174,145],
[251,106,74],
[203,24,29]
]) / 256.

IGNORE_FIRST_PERCENTAGE = 0.1

from math import log10, floor
def round_to_1(x):
	if x < 0:
		x = x*-1
	return -1*round(x, -int(floor(log10(x))))

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	nAvogadro = kb.constants.nAvogadro
	cellDensity = kb.constants.cellDensity

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Calculate concentration data
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	bulkMoleculeIds = bulkMolecules.readAttribute("objectNames")

	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = units.fg * mass.readColumn("cellMass")


	concentrationSetpoints = kb.process.metabolism.metabolitePoolConcentrations
	sortedConcentrationIndex = concentrationSetpoints.asNumber().argsort()[::-1]
	concentrationSetpoints = concentrationSetpoints[sortedConcentrationIndex]

	poolIds = np.array(kb.process.metabolism.metabolitePoolIDs)[sortedConcentrationIndex]
	poolIndexes = np.array([bulkMoleculeIds.index(x) for x in poolIds])
	poolCounts = bulkMolecules.readColumn("counts")[:, poolIndexes]
	poolMols = 1/nAvogadro * poolCounts
	volume = cellMass / cellDensity
	poolConcentrations = poolMols * 1/volume[:,np.newaxis]

	bulkMolecules.close()

	# Compare
	common_units = units.mmol / units.L
	concSetpoint = np.tile(concentrationSetpoints.asNumber(common_units),(time.size,1))
	poolConc = poolConcentrations.asNumber(common_units)
	
	fig = plt.figure(figsize = (11, 11))
	rows = 13
	cols = 12
	idx = 0
	for idx in range(len(poolIds)):
		ax = plt.subplot(rows, cols, idx+1)
		deviation = 1-poolConc[:,idx]/concSetpoint[:,idx]
		ax.plot(time / 60., deviation, linewidth=1, label="pool size", color='k')

		# Highlights >15% deviation
		flag_deviation = abs(deviation[deviation.shape[0] * IGNORE_FIRST_PERCENTAGE:].min())
		bbox = None
		if flag_deviation > 0.15:
			bbox = {'facecolor':'red', 'alpha':0.5, 'pad':1}
		ax.set_title('{}\n{} mmol/L'.format(poolIds[idx][:-3], -1 * round_to_1(concentrationSetpoints[idx].asNumber(units.mmol / units.L))), fontsize=6, bbox=bbox)

		# Sets ticks so that they look pretty
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize=6)
		ax.set_xticks([])
		ymin = deviation[deviation.shape[0] * IGNORE_FIRST_PERCENTAGE:].min()
		ymax = deviation[deviation.shape[0] * IGNORE_FIRST_PERCENTAGE:].max()
		ax.set_ylim([ymin, ymax])
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels([str(round_to_1(deviation.min())), str(round_to_1(deviation.max()))])

	# Create legend
	ax = plt.subplot(rows, cols,len(poolIds) + 1)
	ax.plot(0, 0, linewidth=2, label="1 - c/c_o", color='k')
	ax.legend(loc = 10,prop={'size':10})
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	ax.set_xticks([])
	ax.set_yticks([])
	ax.set_title("Highlights >0.15 deviation", fontsize=12, bbox={'facecolor':'red', 'alpha':0.5, 'pad':1})

	# Save
	plt.subplots_adjust(hspace = 1, wspace = 1)

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
