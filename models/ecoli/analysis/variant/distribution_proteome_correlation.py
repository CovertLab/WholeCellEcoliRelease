from __future__ import absolute_import


import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.utils import parallelization
from wholecell.utils.sparkline import whitePadSparklineAxis
from scipy.stats import pearsonr

SHUFFLE_VARIANT_TAG = "ShuffleParams"
PLACE_HOLDER = -1

FONT_SIZE=9
trim = 0.05


def getPCC((variant, ap, monomerIds, schmidt_counts)):
	try:
		simDir = ap.get_cells(variant = [variant])[0]
		simOutDir = os.path.join(simDir, "simOut")
		sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))

		ids_translation = sim_data.process.translation.monomerData["id"].tolist()
		schmidt_idx = [ids_translation.index(x) for x in monomerIds]

		monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		avgCounts = monomerCounts.readColumn("monomerCounts").mean(axis=0)
		sim_schmidt_counts = avgCounts[schmidt_idx]

		pcc, pval = pearsonr(np.log10(sim_schmidt_counts + 1), np.log10(schmidt_counts + 1))

		return pcc, pval

	except Exception as e:
		print e
		return np.nan, np.nan


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata is not None and SHUFFLE_VARIANT_TAG not in metadata["variant"]:
			print "This plot only runs for variants where parameters are shuffled."
			return

		if not os.path.isdir(inputDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		validation_data = cPickle.load(open(validationDataFile, "rb"))
		schmidtCounts = validation_data.protein.schmidt2015Data["glucoseCounts"]

		ap = AnalysisPaths(inputDir, variant_plot = True)


		pool = parallelization.pool(processes=self.cpus)
		args = zip(range(ap.n_variant), [ap] * ap.n_variant, [validation_data.protein.schmidt2015Data["monomerId"].tolist()] * ap.n_variant, [schmidtCounts] * ap.n_variant)
		result = pool.map(getPCC, args)
		# cPickle.dump(result, open("pcc_results.cPickle", "w"), cPickle.HIGHEST_PROTOCOL)
		pool.close()
		pool.join()
		# result = cPickle.load(open("pcc_results.cPickle", "r"))
		controlPcc, controlPvalue = result[0]
		pccs, pvals = zip(*result[1:])
		pccs = np.array(pccs)
		pvals = np.array(pvals)

		fig = plt.figure()
		fig.set_figwidth(5)
		fig.set_figheight(5)
		ax = plt.subplot(1, 1, 1)

		pccs = np.array([x for x in pccs if not np.isnan(x)])
		ax.hist(pccs, np.sqrt(pccs.size))
		ax.axvline(controlPcc, color = "k", linestyle = "dashed", linewidth = 2)

		ax.set_xlabel("Proteome correlation (Pearson r)")
		ax.set_title("Mean: %0.3g     Std: %0.3g     Control: %0.3g" % (pccs.mean(), pccs.std(), controlPcc))

		axes_list = [ax]

		for a in axes_list:
			for tick in a.yaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)
			for tick in a.xaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)

		whitePadSparklineAxis(ax)

		plt.subplots_adjust(bottom = 0.2, wspace=0.3)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
