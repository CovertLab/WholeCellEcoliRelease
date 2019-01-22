"""
Plot protein monomer counts

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from scipy.stats import pearsonr
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils.fitting import normalize
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get the names of proteins from the KB

		sim_data = cPickle.load(open(simDataFile, "rb"))

		monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
		avgCounts = monomerCounts.readColumn("monomerCounts").mean(axis=0)

		relativeCounts = avgCounts / avgCounts.sum()

		expectedCountsArbitrary = normalize(
			sim_data.process.transcription.rnaExpression[sim_data.condition][sim_data.relation.rnaIndexToMonomerMapping] *
			sim_data.process.translation.translationEfficienciesByMonomer /
			(np.log(2) / sim_data.doubling_time.asNumber(units.s) + sim_data.process.translation.monomerData["degRate"].asNumber(1 / units.s))
			)

		expectedCountsRelative = expectedCountsArbitrary / expectedCountsArbitrary.sum()

		plt.figure(figsize = (8.5, 11))

		maxLine = 1.1 * max(np.log10(expectedCountsRelative.max() + 1), np.log10(relativeCounts.max() + 1))
		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(np.log10(expectedCountsRelative + 1), np.log10(relativeCounts + 1), 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("log10(Expected protein distribution (from fitting))")
		plt.ylabel("log10(Actual protein distribution (average over life cycle))")
		plt.title("PCC (of log values): %0.2f" % pearsonr(np.log10(expectedCountsRelative + 1), np.log10(relativeCounts + 1))[0])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
