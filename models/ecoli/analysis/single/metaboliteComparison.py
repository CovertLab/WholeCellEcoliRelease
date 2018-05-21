#!/usr/bin/env python
"""
Compare metabolite concentrations to target concentrations in FBA objective

@date: Created 2/6/18
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils.sparkline import whitePadSparklineAxis


def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, 'simOutDir does not currently exist as a directory'

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Read data from listeners
	enzymeKinetics = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	actualConc = enzymeKinetics.readColumn("metaboliteConcentrations")[1:,:]
	enzymeKinetics.close()

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	targetConc = fbaResults.readColumn("targetConcentrations")[1:,:]
	moleculeNames = np.array(fbaResults.readAttribute("homeostaticTargetMolecules"))
	fbaResults.close()

	mainReader = TableReader(os.path.join(simOutDir, "Main"))
	time = mainReader.readColumn("time")[1:] - mainReader.readAttribute("initialTime")
	mainReader.close()

	# Average concentrations over time for subplot 1
	actualConcAve = np.nanmean(actualConc, axis=0)
	targetConcAve = np.nanmean(targetConc, axis=0)

	plt.figure(figsize = (8.5, 11))

	# Plot average comparison with lines denoting order of magnitude
	ax = plt.subplot(3,1,1)
	ax.plot([-6, 6], [-6, 6], 'k')
	ax.plot([-5, 6], [-6, 5], 'k')
	ax.plot([-6, 5], [-5, 6], 'k')
	ax.plot(np.log10(targetConcAve), np.log10(actualConcAve), "ob", markeredgewidth=0, alpha=0.25)
	ax.set_xlabel("Log10(Target Concentration [mmol/L])", fontsize=8)
	ax.set_ylabel("Log10(Actual Concentration [mmol/L])", fontsize=8)

	whitePadSparklineAxis(ax)
	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	ax.set_ylim(ylim[0] - 0.1, ylim[1])
	ax.set_xlim(xlim[0] - 0.1, xlim[1])
	ax.set_yticks(range(-6, int(ylim[1]) + 1, 2))
	ax.set_xticks(range(-6, int(xlim[1]) + 1, 2))
	ax.tick_params(axis='both', which='major', labelsize=6)

	# Plot ratio of actual concentration to target concentration
	ratio = np.log10(actualConc / targetConc)
	ax = plt.subplot(3,1,2)
	ax.plot(time / 60, ratio)
	ax.set_xlabel("Time (min)", fontsize=8)
	ax.set_ylabel("Log10(Concentration to Target)", fontsize=8)
	ax.tick_params(axis='both', which='major', labelsize=6)

	# Plot outliers of ratio
	means = np.mean(ratio, axis=0)
	mean = np.mean(means)
	std = np.std(means)
	outliers = np.unique(np.where((ratio[1:, :] > mean + std / 2) | (ratio[1:, :] < mean - std / 2))[1])
	idx = outliers[np.argsort(means[outliers])][::-1]

	ax = plt.subplot(3,1,3)
	if len(idx):
		ax.plot(time / 60, ratio[:, idx])
		ax.legend(moleculeNames[idx], fontsize=6)
	ax.set_xlabel("Time (min)", fontsize=8)
	ax.set_ylabel("Log10(Concentration to Target)", fontsize=8)
	ax.tick_params(axis='both', which='major', labelsize=6)

	plt.tight_layout()

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")


if __name__ == '__main__':
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument('simOutDir', help = 'Directory containing simulation output', type = str)
	parser.add_argument('plotOutDir', help = 'Directory containing plot output (will get created if necessary)', type = str)
	parser.add_argument('plotOutFileName', help = 'File name to produce', type = str)
	parser.add_argument('--simDataFile', help = 'KB file name', type = str, default = defaultSimDataFile)
	parser.add_argument('--validationDataFile')

	args = parser.parse_args().__dict__
	main(args['simOutDir'], args['plotOutDir'], args['plotOutFileName'], args['simDataFile'], args['validationDataFile'])
