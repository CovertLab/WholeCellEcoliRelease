"""
Plot first-order rate constants of mRNAs, observed vs expected.
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

# Observed degradation rates are only calculated for RNAs with mean counts
# no less than this threshold
MEAN_RNA_COUNT_THRESHOLD = 3


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get the expected degradation rates from KB
		sim_data = cPickle.load(open(simDataFile, 'rb'))
		mRNA_ids = sim_data.process.transcription.rna_data['id']
		isMRna = sim_data.process.transcription.rna_data['is_mRNA']
		expected_degradation_rate_constants = np.array(
			sim_data.process.transcription.rna_data['deg_rate'][isMRna].asNumber()
			)

		# Get length of simulation
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		sim_length = main_reader.readColumn('time')[-1]

		# Read counts of mRNAs
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		mRNA_counts_mean = mRNA_counts.mean(axis=0)

		# Read number of degradation events
		rna_degradation_reader = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		n_RNA_degraded = rna_degradation_reader.readColumn('countRnaDegraded')
		n_mRNA_degraded_total = n_RNA_degraded.sum(axis = 0)[isMRna]

		# Get mask for RNAs with counts no less than threshold
		mask = mRNA_counts_mean >= MEAN_RNA_COUNT_THRESHOLD

		if mask.sum() == 0:
			print("Skipping analysis - RNA counts not sufficient for analysis")
			return

		# Calculate observed rate constants
		observed_rate_constants = (n_mRNA_degraded_total[mask] / sim_length) / mRNA_counts_mean[mask]
		expected_rate_constants = expected_degradation_rate_constants[mask]

		plt.figure(figsize = (8, 8))
		maxLine = 1.1*np.max(
			np.concatenate((observed_rate_constants, expected_rate_constants)))

		plt.plot([0, maxLine], [0, maxLine], '--r')
		plt.plot(expected_rate_constants, observed_rate_constants,
			'o', markeredgecolor = 'k', markerfacecolor = 'none')

		plt.xlabel("Expected RNA decay rate constant [$s^{-1}$]")
		plt.ylabel("Observed RNA decay rate constant [$s^{-1}$]")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
