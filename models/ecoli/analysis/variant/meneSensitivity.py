"""
Plot to assess sensitivity of menE behavior to model parameters.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/1/17
"""

from __future__ import absolute_import


import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

THRESHOLD = 0.001 # .1 percent
TARGET_CONC_SINGLE = 0.10183094010881857 * units.mmol / units.L # found from WT sim; mmol/L
FONTSIZE = 10
LABELSIZE = 8
MARKERSIZE = 1


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, "inputDir does not currently exist as a directory"
		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get cells
		ap = AnalysisPaths(inputDir, variant_plot = True)
		if ap.n_variant != 9:
			print "This plot expects all variants of mene_params"
			return

		# Get constants from wildtype variant
		sim_data = cPickle.load(open(ap.get_variant_kb(4), "rb")) # 4 is the wildtype variant
		cellDensity = sim_data.constants.cellDensity
		nAvogadro = sim_data.constants.nAvogadro

		# Initialize variables
		enzymeId = "MENE-CPLX[c]"
		endProductIds = ["REDUCED-MENAQUINONE[c]", "CPD-12115[c]"]
		TARGET_CONC = len(endProductIds) * TARGET_CONC_SINGLE

		# Investigate each variant
		meneDepletion = np.zeros([ap.n_seed, ap.n_variant])
		endProductDepletion = np.zeros([ap.n_seed, ap.n_variant])

		simOutDir = None

		for variant in xrange(ap.n_variant):
			for seed in xrange(ap.n_seed):
				cells = ap.get_cells(variant = [variant], seed = [seed])
				timeMeneDepleted = [] # seconds
				timeEndProdDepleted = [] # seconds

				for i, simDir in enumerate(cells):
					simOutDir = os.path.join(simDir, "simOut")

					# Get molecule counts
					bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
					moleculeIds = bulkMolecules.readAttribute("objectNames")
					bulkMoleculeCounts = bulkMolecules.readColumn("counts")

					meneIndex = moleculeIds.index(enzymeId)
					meneCounts = bulkMoleculeCounts[:, meneIndex]
					endProductIndices = [moleculeIds.index(x) for x in endProductIds]
					endProductCounts = bulkMoleculeCounts[:, endProductIndices]
					bulkMolecules.close()

					# Compute time with zero counts of tetramer (MENE-CPLX)
					timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
					meneDepletionIndices = np.where(meneCounts == 0)[0]
					timeMeneDepleted.append(timeStepSec[meneDepletionIndices].sum())

					# Compute time with end products under the target concentration
					mass = TableReader(os.path.join(simOutDir, "Mass")).readColumn("cellMass") * units.fg
					volume = mass / cellDensity
					endProductConcentrations = np.sum([endProductCounts[:, col] / nAvogadro / volume for col in xrange(endProductCounts.shape[1])], axis = 0)
					endProductDepletionIndices = np.where(endProductConcentrations < ((1 - THRESHOLD) * TARGET_CONC))[0]
					timeEndProdDepleted.append(timeStepSec[endProductDepletionIndices].sum())

				# Record MENE-CPLX depletion
				# TODO(jerry): Is this meant to read from the inner loop's last simOutDir?
				totalTime = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")[-1] + timeStepSec[-1]
				fractionMeneDepleted = np.sum(timeMeneDepleted) / totalTime
				meneDepletion[seed, variant] = fractionMeneDepleted

				# Record end product depletion
				fractionEndProdDepleted = np.sum(timeEndProdDepleted) / totalTime
				endProductDepletion[seed, variant] = fractionEndProdDepleted

		# Compute average and standard deviations
		meneDepletion_avg = np.average(meneDepletion, axis = 0)
		meneDepletion_std = np.std(meneDepletion, axis = 0)
		endProductDepletion_avg = np.average(endProductDepletion, axis = 0)
		endProductDepletion_std = np.std(endProductDepletion, axis = 0)

		# Plot
		fig, axesList = plt.subplots(2, 1, figsize = (8, 8))
		ax1, ax2 = axesList
		xvals = np.arange(ap.n_variant)
		fig.suptitle("Sensitivity Analysis: menE depletion")

		for ax, avg, std in zip(axesList, [meneDepletion_avg, endProductDepletion_avg], [meneDepletion_std, endProductDepletion_std]):
			ax.scatter(xvals, avg, edgecolor = "none", clip_on = False, s = MARKERSIZE)
			ax.errorbar(xvals, avg, yerr = std, color = "b", linewidth = 1, clip_on = False, fmt = "o", capsize = 4, capthick = 1, markeredgecolor = "none")

		ax1.set_title("MenE tetramer depletion", fontsize = FONTSIZE)
		ax2.set_title("Menaquinone products depletion", fontsize = FONTSIZE)
		xlabels = ["1/10 x", "1/8 x", "1/4 x", "1/2 x", "1 x", "2 x", "4 x", "8 x", "10 x"]
		title_tags = ["counts = 0", "<%s percent of wildtype" % (THRESHOLD * 100)]
		for i, ax in enumerate([ax1, ax2]):
			ax.set_ylabel("Fraction of Time\n%s" % title_tags[i], fontsize = FONTSIZE)
			ax.set_xlabel("Factor of increase of menE synthesis probability", fontsize = FONTSIZE)
			ax.set_ylim([-0.1, 1.1])
			ax.set_xlim([-0.25, 8.25])
			whitePadSparklineAxis(ax)
			ax.set_xticks(xvals)
			ax.set_xticklabels(xlabels)
			ax.set_yticks([0, 1])

		plt.subplots_adjust(hspace = 1, wspace = 1, top = 0.9, bottom = 0.1)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		# Plot clean versions for figure
		FIRST = True
		for avg, std, filename in zip([meneDepletion_avg, endProductDepletion_avg], [meneDepletion_std, endProductDepletion_std], ["mene", "menaquinone"]):
			fig, ax = plt.subplots(1, 1, figsize = (10, 3))
			ax.scatter(xvals, avg, edgecolor = "none", clip_on = False, s = MARKERSIZE)
			ax.errorbar(xvals, avg, yerr = std, color = "b", linewidth = 1, clip_on = False, fmt = "o", capsize = 4, capthick = 1, markeredgecolor = "none")
			ax.set_ylim([-0.1, 1.1])
			ax.set_xlim([-0.25, 8.25])
			if FIRST:
				FIRST = False
				whitePadSparklineAxis(ax, False)
			else:
				whitePadSparklineAxis(ax)
			ax.set_xticks(xvals)
			ax.set_xticklabels([])
			ax.set_yticks([0, 1])
			ax.set_yticklabels([])
			exportFigure(plt, plotOutDir, plotOutFileName + "_%s" % filename, metadata)
			plt.close("all")


if __name__ == "__main__":
	Plot().cli()
