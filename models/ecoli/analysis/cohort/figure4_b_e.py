#!/usr/bin/env python

import argparse
import os
import re

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

from wholecell.utils.sparkline import whitePadSparklineAxis


PLACE_HOLDER = -1

FONT_SIZE=8
trim = 0.03


# def sparklineAxis(axis):
# 	axis.spines['top'].set_visible(False)
# 	axis.spines['bottom'].set_visible(False)
# 	axis.xaxis.set_ticks_position('none')
# 	axis.tick_params(which = 'both', direction = 'out')


def mm2inch(value):
	return value * 0.0393701

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile = None, metadata = None):

	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	import cPickle
	sim_data = cPickle.load(open(simDataFile, "rb"))
	oriC = sim_data.constants.oriCCenter.asNumber()
	terC = sim_data.constants.terCCenter.asNumber()
	genomeLength = len(sim_data.process.replication.genome_sequence)


	mult = 2.3
	fig = plt.figure()
	fig.set_figwidth(mm2inch(97)*mult)
	fig.set_figheight(mm2inch(58)*mult)

	ax0 = plt.subplot2grid((3,4), (0,0), colspan = 4)
	ax1 = plt.subplot2grid((3,4), (1,0), colspan = 4, sharex=ax0)
	ax2 = plt.subplot2grid((3,4), (2,0), colspan = 4, sharex=ax0)
	ax0_xlim = ax0.get_xlim()
	#ax3 = plt.subplot2grid((5,7), (3,0), colspan = 4, sharex=ax0)
	#ax4 = plt.subplot2grid((5,7), (4,0), colspan = 4, sharex=ax0)

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	#all_cells = ap.get_cells(seed=[0,1,2,3])
	all_cells = ap.get_cells(seed=[4])

	if not len(all_cells):
		return

	for idx, simDir in enumerate(all_cells):
		color = "black"
		alpha = 1
		if idx % 2:
			color = "#BF673B"
			blue = 0.8

		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		## Cell mass
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = mass.readColumn("cellMass")
		massPerOric = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")
		idxInit = np.where(massPerOric >= 1)[0]
		ax0.plot(time / 60., cellMass, color = color, alpha = alpha, linewidth=2)
		ax0.plot(time[idxInit] / 60., cellMass[idxInit],  markersize=6, linewidth=0, marker="o", color = "#FCBE67", markeredgewidth=0)

		## Inst. growth rate
		growthRate = mass.readColumn("instantaniousGrowthRate")
		growthRate = (1 / units.s) * growthRate
		growthRate = growthRate.asNumber(1 / units.min)
		growthRate[abs(growthRate - np.median(growthRate)) > 1.25 * np.nanstd(growthRate)] = np.nan
		ax1.plot(time / 60., growthRate, color = color, alpha = alpha)

		## Rna over protein
		# Get active ribosome counts
		uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		ribosomeCounts = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
		uniqueMoleculeCounts.close()
		ribosomeConcentration = ((1 / sim_data.constants.nAvogadro) * ribosomeCounts) / ((1.0 / sim_data.constants.cellDensity) * (units.fg * cellMass))
		ribosomeConcentration = ribosomeConcentration.asNumber(units.umol / units.L)
		ax2.plot(time / 60., ribosomeConcentration, color = color, alpha = alpha, linewidth=2)
		ax2.set_ylim([18., 28.])
		# rnaMass = mass.readColumn("rnaMass")
		# proteinMass = mass.readColumn("proteinMass")
		# rnaOverProtein = rnaMass / proteinMass
		# ax2.plot(time / 60., rnaOverProtein, color = color, alpha = alpha, linewidth=2)


		## Fork position
		sequenceIdx = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceIdx")
		sequenceLength = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceLength")
		reverseIdx = 1
		reverseCompIdx = 3
		reverseSequences = np.logical_or(sequenceIdx == reverseIdx, sequenceIdx == reverseCompIdx)
		sequenceLength[reverseSequences] = -1 * sequenceLength[reverseSequences]
		sequenceLength[sequenceLength == PLACE_HOLDER] = np.nan

		# Down sample dna polymerase position, every position is only plotted once here
		# using numpy ninja-ness
		unique, index, value = np.unique(sequenceLength, return_index=True, return_inverse=True)
		m = np.zeros_like(value, dtype=bool)
		m[index] = True
		m = m.reshape(sequenceLength.shape)
		sequenceLength[~m] = np.nan

		# ax3.plot(time / 60., sequenceLength, marker=',', markersize=2, linewidth=0, color = color, alpha = alpha)
		# ax3.set_yticks([-1 * genomeLength / 2, 0, genomeLength / 2])
		# ax3.set_yticklabels(['-terC', 'oriC', '+terC'])

		## Pairs of forks
		# pairsOfForks = (sequenceIdx != PLACE_HOLDER).sum(axis = 1) / 4
		# ax4.plot(time / 60., pairsOfForks, linewidth=2, color = color, alpha = alpha)
		# ax4.set_yticks(np.arange(0,7))
		# ax4.set_ylim([0, 6])


	y_ticks = ax0.get_yticks()
	new_y_ticks = y_ticks[0:-1:2]
	ax0.set_yticks(new_y_ticks)

	ax0.set_xlim([0., time.max() / 60.])
	ax0.set_ylabel("Cell mass\n(fg)", fontsize=FONT_SIZE)
	ax0.xaxis.set_visible(False)
	#ax0.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
	try:
		T_ADD_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][1][0] / 60.
	except:
		print "nutrientsTimeSeries does not have correct dimensions for this analysis. Exiting."
		return
	axes_list = [ax0, ax1, ax2]#, ax3, ax4]
	for a in axes_list:
		shift_time = T_ADD_AA
		width = a.get_xlim()[1] - shift_time
		height = a.get_ylim()[1] - a.get_ylim()[0]
		a.add_patch(
			patches.Rectangle(
				(shift_time, a.get_ylim()[0]),   # (x,y)
				width,          # width
				height,          # height
				alpha = 0.25,
				color = "gray",
				linewidth = 0.
			)
		)

	for a in axes_list:
		for tick in a.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

	ax1.set_ylabel(r"$\mu$ $(\frac{gDCW}{gDCW \cdot \, min})$", fontsize=FONT_SIZE)
	ax1.xaxis.set_visible(False)
	# ax1.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)
	ax1.set_ylim([0.008, 0.032])

	#ax2.set_ylabel("RNA/Protein\n(fg/fg)", fontsize=FONT_SIZE)
	ax2.set_ylabel("Active\nribosome\n(umol/L)", fontsize=FONT_SIZE)
	ax2.xaxis.set_visible(False)
	# ax2.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	y_ticks = ax2.get_yticks()
	new_y_ticks = y_ticks[0:-1:2]
	ax2.set_yticks(new_y_ticks)

	# ax3.set_ylabel("DNA polymerase\nposition (nt)", fontsize=FONT_SIZE)
	# ax3.xaxis.set_visible(False)
	# # ax3.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	# ax4.set_ylabel("Relative rate\nof dNTP\npolymerization", fontsize=FONT_SIZE)
	# ax4.set_xlabel("Time (min)", fontsize=FONT_SIZE)
	# # ax4.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	whitePadSparklineAxis(ax0, False)
	whitePadSparklineAxis(ax1, False)
	whitePadSparklineAxis(ax2, False)
	# whitePadSparklineAxis(ax3, False)
	# whitePadSparklineAxis(ax4)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_b", metadata)


	for axes in axes_list:
		axes.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom='off',      # ticks along the bottom edge are off
			top='off',         # ticks along the top edge are off
			labelbottom='off') # labels along the bottom edge are off
		axes.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			left='off',      # ticks along the bottom edge are off
			right='off',         # ticks along the top edge are off
			labelleft='off') # labels along the bottom edge are off

		axes.set_xlabel("")
		axes.set_ylabel("")

	plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_b_stripped" ,metadata, transparent = True)
	plt.close("all")

	################ PART II #######################

	mult = 2.3
	fig = plt.figure()
	fig.set_figwidth(mm2inch(97)*mult)
	fig.set_figheight(mm2inch(38)*mult)

	# ax0 = plt.subplot2grid((3,4), (0,0), colspan = 4)
	# ax1 = plt.subplot2grid((3,4), (1,0), colspan = 4, sharex=ax0)
	# ax2 = plt.subplot2grid((3,4), (2,0), colspan = 4, sharex=ax0)
	ax3 = plt.subplot2grid((2,4), (0,0), colspan = 4, sharex = ax0)
	ax4 = plt.subplot2grid((2,4), (1,0), colspan = 4, sharex=ax0)

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	all_cells = ap.get_cells(seed=[4])

	if not len(all_cells):
		return

	for idx, simDir in enumerate(all_cells):
		color = "black"
		alpha = 1
		if idx % 2:
			color = "#BF673B"
			blue = 0.8

		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

		## Cell mass
		# mass = TableReader(os.path.join(simOutDir, "Mass"))
		# cellMass = mass.readColumn("cellMass")
		# ax0.plot(time / 60., cellMass, color = color, alpha = alpha, linewidth=2)

		## Inst. growth rate
		# growthRate = mass.readColumn("instantaniousGrowthRate")
		# growthRate = (1 / units.s) * growthRate
		# growthRate = growthRate.asNumber(1 / units.min)
		# growthRate[abs(growthRate - np.median(growthRate)) > 1.25 * np.nanstd(growthRate)] = np.nan
		# ax1.plot(time / 60., growthRate, color = color, alpha = alpha)

		## Rna over protein
		# Get active ribosome counts
		# uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		# ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
		# ribosomeCounts = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]
		# uniqueMoleculeCounts.close()
		# ribosomeConcentration = ((1 / sim_data.constants.nAvogadro) * ribosomeCounts) / ((1.0 / sim_data.constants.cellDensity) * (units.fg * cellMass))
		# ribosomeConcentration = ribosomeConcentration.asNumber(units.umol / units.L)
		# ax2.plot(time / 60., ribosomeConcentration, color = color, alpha = alpha, linewidth=2)
		# ax2.set_ylim([20., 35.])
		# rnaMass = mass.readColumn("rnaMass")
		# proteinMass = mass.readColumn("proteinMass")
		# rnaOverProtein = rnaMass / proteinMass
		# ax2.plot(time / 60., rnaOverProtein, color = color, alpha = alpha, linewidth=2)


		## Fork position
		sequenceIdx = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceIdx")
		sequenceLength = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceLength")
		reverseIdx = 1
		reverseCompIdx = 3
		reverseSequences = np.logical_or(sequenceIdx == reverseIdx, sequenceIdx == reverseCompIdx)
		sequenceLength[reverseSequences] = -1 * sequenceLength[reverseSequences]
		sequenceLength[sequenceLength == PLACE_HOLDER] = np.nan

		# Down sample dna polymerase position, every position is only plotted once here
		# using numpy ninja-ness
		unique, index, value = np.unique(sequenceLength, return_index=True, return_inverse=True)
		m = np.zeros_like(value, dtype=bool)
		m[index] = True
		m = m.reshape(sequenceLength.shape)
		sequenceLength[~m] = np.nan

		ax3.plot(time / 60., sequenceLength, marker=',', markersize=2, linewidth=0, color = color, alpha = alpha)
		ax3.set_yticks([-1 * genomeLength / 2, 0, genomeLength / 2])
		ax3.set_yticklabels(['-terC', 'oriC', '+terC'])

		# Pairs of forks
		pairsOfForks = (sequenceIdx != PLACE_HOLDER).sum(axis = 1) / 4
		ax4.plot(time / 60., pairsOfForks, linewidth=2, color = color, alpha = alpha)
		ax4.set_yticks(np.arange(0,7))
		ax4.set_ylim([0, 6.1])
		ax4.set_yticks([0, 6.])
		ax4.set_yticklabels(["0","6"])


	# y_ticks = ax0.get_yticks()
	# new_y_ticks = y_ticks[0:-1:2]
	# ax0.set_yticks(new_y_ticks)

	ax3.set_xlim([0., time.max() / 60.])
	# ax0.set_ylabel("Cell mass\n(fg)", fontsize=FONT_SIZE)
	# ax0.xaxis.set_visible(False)
	#ax0.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	nutrientsTimeSeriesLabel = sim_data.nutrientsTimeSeriesLabel
	T_ADD_AA = sim_data.nutrientsTimeSeries[nutrientsTimeSeriesLabel][1][0] / 60.
	axes_list = [ax3, ax4]
	for a in axes_list:
		shift_time = T_ADD_AA
		width = a.get_xlim()[1] - shift_time
		height = a.get_ylim()[1] - a.get_ylim()[0]
		a.add_patch(
			patches.Rectangle(
				(shift_time, a.get_ylim()[0]),   # (x,y)
				width,          # width
				height,          # height
				alpha = 0.25,
				color = "gray",
				linewidth = 0.
			)
		)

	for a in axes_list:
		for tick in a.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

	# ax1.set_ylabel(r"$\mu$ $(\frac{gDCW}{gDCW-min})$", fontsize=FONT_SIZE)
	# ax1.xaxis.set_visible(False)
	# # ax1.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)
	# ax1.set_ylim([0.012, 0.032])

	# #ax2.set_ylabel("RNA/Protein\n(fg/fg)", fontsize=FONT_SIZE)
	# ax2.set_ylabel("Active\nribosome\n(umol/L)", fontsize=FONT_SIZE)
	# ax2.xaxis.set_visible(False)
	# # ax2.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	# y_ticks = ax2.get_yticks()
	# new_y_ticks = y_ticks[0:-1:2]
	# ax2.set_yticks(new_y_ticks)

	ax3.set_ylabel("DNA polymerase\nposition (nt)", fontsize=FONT_SIZE)
	ax3.xaxis.set_visible(False)
	# ax3.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	ax4.set_ylabel("Relative rate\nof dNTP\npolymerization", fontsize=FONT_SIZE)
	ax4.set_xlabel("Time (min)", fontsize=FONT_SIZE)
	# ax4.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	# whitePadSparklineAxis(ax0, False)
	# whitePadSparklineAxis(ax1, False)
	# whitePadSparklineAxis(ax2, False)
	whitePadSparklineAxis(ax3, False)
	whitePadSparklineAxis(ax4)

	plt.subplots_adjust(bottom=0.15)


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_e", metadata)


	for axes in axes_list:
		axes.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom='off',      # ticks along the bottom edge are off
			top='off',         # ticks along the top edge are off
			labelbottom='off') # labels along the bottom edge are off
		axes.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			left='off',      # ticks along the bottom edge are off
			right='off',         # ticks along the top edge are off
			labelleft='off') # labels along the bottom edge are off

		axes.set_xlabel("")
		axes.set_ylabel("")

	plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_e_stripped" ,metadata, transparent = True)
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
