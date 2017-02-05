#!/usr/bin/env python

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

PLACE_HOLDER = -1

FONT_SIZE=9

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


	fig = plt.figure()
	fig.set_figwidth(8.5)
	fig.set_figheight(5.5)

	ax0 = plt.subplot2grid((5,7), (0,0), colspan = 4)
	ax1 = plt.subplot2grid((5,7), (1,0), colspan = 4, sharex=ax0)
	ax2 = plt.subplot2grid((5,7), (2,0), colspan = 4, sharex=ax0)
	ax3 = plt.subplot2grid((5,7), (3,0), colspan = 4, sharex=ax0)
	ax4 = plt.subplot2grid((5,7), (4,0), colspan = 4, sharex=ax0)

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)
	all_cells = ap.get_cells(seed=[0])

	for idx, simDir in enumerate(all_cells):
		color = "black"
		alpha = 0.8
		if idx % 2:
			color = "blue"
			blue = 0.8

		simOutDir = os.path.join(simDir, "simOut")

		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
		
		## Cell mass
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = mass.readColumn("cellMass")
		ax0.plot(time / 60., cellMass, color = color, alpha = alpha, linewidth=2)

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
		ax2.set_ylim([20., 35.])
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

		## Pairs of forks
		pairsOfForks = (sequenceIdx != PLACE_HOLDER).sum(axis = 1) / 4
		ax4.plot(time / 60., pairsOfForks, linewidth=2, color = color, alpha = alpha)
		ax4.set_yticks(np.arange(0,7))
		ax4.set_ylim([0, 6])


	ax0.set_xlim([0., time.max() / 60.])
	ax0.set_ylabel("Cell mass\n(fg)", fontsize=FONT_SIZE)
	ax0.xaxis.set_visible(False)
	ax0.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	ax1.set_ylabel(r"$\mu$ $(\frac{gDCW}{gDCW-min})$", fontsize=FONT_SIZE)
	ax1.xaxis.set_visible(False)
	ax1.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)
	ax1.set_ylim([0.012, 0.032])

	#ax2.set_ylabel("RNA/Protein\n(fg/fg)", fontsize=FONT_SIZE)
	ax2.set_ylabel("Active\nribosome\n(umol/L)", fontsize=FONT_SIZE)
	ax2.xaxis.set_visible(False)
	ax2.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	ax3.set_ylabel("DNA polymerase\nposition (nt)", fontsize=FONT_SIZE)
	ax3.xaxis.set_visible(False)
	ax3.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)

	ax4.set_ylabel("Relative rate\nof dNTP\npolymerization", fontsize=FONT_SIZE)
	ax4.xaxis.set_visible(False)
	ax4.axvline(x=44*2+22., linewidth=3, color='gray', alpha = 0.5)


	ax5 = plt.subplot2grid((5,7), (0,5), colspan = 2, rowspan = 2)

	ax5.set_ylabel("Added mass (fg)")
	ax5.set_xlabel("Initial mass (fg)")


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
