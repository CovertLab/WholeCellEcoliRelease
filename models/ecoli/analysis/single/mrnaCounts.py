"""
Plot mRNA counts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Get the names of RNAs from the KB
		is_mRNA = sim_data.process.transcription.rnaData['isMRna']
		mRNA_ids = sim_data.process.transcription.rnaData['id'][is_mRNA]

		# Get reader for mRNA counts
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))

		# Check that the order of mRNAs in table matches that of KB
		assert np.all(mRNA_ids == mRNA_counts_reader.readAttribute('mRNA_ids'))

		# Read final mRNA counts from reader
		counts = mRNA_counts_reader.readColumn('mRNA_counts')[-1, :]

		plt.figure(figsize = (8.5, 11))

		expectedCountsArbitrary = sim_data.process.transcription.rnaExpression[
			sim_data.condition][is_mRNA]
		expectedCounts = expectedCountsArbitrary/expectedCountsArbitrary.sum() * counts.sum()

		maxLine = 1.1 * max(expectedCounts.max(), counts.max())
		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(expectedCounts, counts, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("Expected RNA count (scaled to total)")
		plt.ylabel("Actual RNA count (at final time step)")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
