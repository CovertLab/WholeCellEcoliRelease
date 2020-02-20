"""
Plots various effects that may be limiting growth

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/2015
"""

from __future__ import absolute_import, division, print_function

import os

from math import log10, floor
import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot

COLOR_CHOICES = np.array([
	[0,0,0],
	[252,174,145],
	[251,106,74],
	[203,24,29]
]) / 256.

IGNORE_FIRST_PERCENTAGE = 0.1

TOP_RANGE_MARKER_COLOR = 'red'
BOTTOM_RANGE_MARKER_COLOR = 'blue'
BOTH_RANGE_MARKER_COLOR = 'purple'


def round_to_1(x):
	if x < 0:
		x = x*-1
	return -1*round(x, -int(floor(log10(x))))

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		cellDensity = sim_data.constants.cellDensity
		homeostaticRangeObjFractionHigher = sim_data.constants.metabolismHomeostaticRangeObjFractionHigher

		# Load time
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Calculate concentration data
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = units.fg * mass.readColumn("cellMass")

		concIds = sorted(sim_data.process.metabolism.concDict)
		concPools = units.mol / units.L * np.array([sim_data.process.metabolism.concDict[key].asNumber(units.mol / units.L) for key in concIds])
		concentrationSetpoints = concPools
		sortedConcentrationIndex = concentrationSetpoints.asNumber().argsort()[::-1]
		concentrationSetpoints = concentrationSetpoints[sortedConcentrationIndex]

		poolIds = np.array(concIds)[sortedConcentrationIndex]
		(poolCounts,) = read_bulk_molecule_counts(simOutDir, (poolIds,))
		poolMols = 1/nAvogadro * poolCounts
		volume = cellMass / cellDensity
		poolConcentrations = poolMols * 1/volume[:,np.newaxis]

		# Compare
		common_units = units.mmol / units.L
		concSetpoint = np.tile(concentrationSetpoints.asNumber(common_units),(time.size,1))
		poolConc = poolConcentrations.asNumber(common_units)

		plt.figure(figsize=(15, 15))
		n_mets = len(poolIds)
		n_subplots = n_mets + 2  # need additional subplots for the legend
		rows = int(np.ceil(np.sqrt(n_subplots)))
		cols = int(np.ceil(n_subplots / rows))
		for idx in range(n_mets):
			ax = plt.subplot(rows, cols, idx+1)
			deviation = 1-poolConc[:,idx]/concSetpoint[:,idx]
			ax.plot(time / 60., deviation, linewidth=1, label="pool size", color='k')

			# Highlights >15% deviation
			flag_deviation = np.abs(deviation[int(deviation.shape[0] * IGNORE_FIRST_PERCENTAGE):]).max()
			bbox = None
			if flag_deviation > 0.15:
				bbox = {'facecolor':'red', 'alpha':0.5, 'pad':1}
			ax.set_title('{}\n{:.2g} mmol/L'.format(poolIds[idx][:-3], concentrationSetpoints[idx].asNumber(units.mmol / units.L)), fontsize=6, bbox=bbox)

			# Sets ticks so that they look pretty
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize=6)
			ax.set_xticks([])

			targetRange = (0, homeostaticRangeObjFractionHigher)
			targetMin = np.amin(targetRange)
			targetMax = np.amax(targetRange)

			# Put lines for above and below target range
			if targetMin == targetMax:
				ax.axhline(y=targetMax, color=BOTH_RANGE_MARKER_COLOR)
			else:
				ax.axhline(y=targetMax, color=TOP_RANGE_MARKER_COLOR)
				ax.axhline(y=targetMin, color=BOTTOM_RANGE_MARKER_COLOR)
			ymin = deviation[int(deviation.shape[0] * IGNORE_FIRST_PERCENTAGE):].min()
			ymax = deviation[int(deviation.shape[0] * IGNORE_FIRST_PERCENTAGE):].max()
			if ymax < 0:
				ymax = 0
			if ymin > 0:
				ymin = 0

			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(['{:.2g}'.format(ymin), '{:.2g}'.format(ymax)])

		# Create legend
		ax = plt.subplot(rows, cols, n_subplots)
		ax.plot(0, 0, linewidth=1, label="1 - c/c_o", color='k')
		if homeostaticRangeObjFractionHigher == 0:
			ax.plot(0, 0, linewidth=1, label="Target", color=BOTH_RANGE_MARKER_COLOR)
		else:
			ax.plot(0, 0, linewidth=1, label="Top of target range", color=TOP_RANGE_MARKER_COLOR)
			ax.plot(0, 0, linewidth=1, label="Bottom of target range", color=BOTTOM_RANGE_MARKER_COLOR)
		ax.legend(loc = 'lower center',prop={'size':"xx-small"})
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.yaxis.set_ticks_position('none')
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_title("Highlights >0.15 deviation", fontsize='x-small', bbox={'facecolor':'red', 'alpha':0.5, 'pad':1})

		# Save
		plt.subplots_adjust(hspace = 1, wspace = 1)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
