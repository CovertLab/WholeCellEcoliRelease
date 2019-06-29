"""
Computes the average distance between transcribing RNA polymerases for each
gene (transcription unit), and compares the distance to the known size of the
molecular footprint of RNAP.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/25/19
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath, units


SAMPLE_SIZE = 75  # Number of genes to plot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		rnap_data_reader = TableReader(os.path.join(simOutDir, 'RnapData'))
		rna_synth_prob_reader = TableReader(os.path.join(simOutDir, 'RnaSynthProb'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = (units.s)*(main_reader.readColumn('time') - initial_time)
		n_rna_init_events = rnap_data_reader.readColumn('rnaInitEvent')
		gene_copy_numbers = rna_synth_prob_reader.readColumn('gene_copy_number')
		TU_ids = rna_synth_prob_reader.readAttribute('rnaIds')

		# Get RNAP elongation rate for sim condition
		nutrients = sim_data.conditions[sim_data.condition]["nutrients"]
		rnap_elong_rate = sim_data.process.transcription.rnaPolymeraseElongationRateDict[nutrients]

		# Calculate the total number of initiation events that happen to each
		# gene per gene copy throughout the cell cycle
		n_total_rna_init_events_per_copy = np.divide(
			n_rna_init_events[1:, :].astype(np.float64),
			gene_copy_numbers[1:, :]
			).sum(axis=0)

		# Divide by length of cell cycle to get average initiation rate
		avg_init_rate = (1./time[-1])*n_total_rna_init_events_per_copy

		# Divide elongation rate with initiation rate to get average distance
		# between RNAPs in nucleotides
		avg_inter_rnap_distance = (rnap_elong_rate/avg_init_rate).asNumber(units.nt)

		# Sort from shortest to longest
		sorted_index = avg_inter_rnap_distance.argsort()
		sorted_TU_ids = [TU_ids[i][:-3] for i in sorted_index]  # [c] stripped
		avg_inter_rnap_distance.sort()

		# Get RNAP footprint size from validation data
		RNAP_footprint_size = validation_data.dna_footprint_sizes[
			sim_data.moleculeIds.rnapFull].asNumber(units.nt)

		# Mark genes with RNAPs that are too close to each other
		n_too_close = (avg_inter_rnap_distance[:SAMPLE_SIZE] < RNAP_footprint_size).sum()
		bar_colors = ["r"]*n_too_close + ["b"]*(SAMPLE_SIZE - n_too_close)

		# Plot the first n genes with shortest distances
		plt.figure(figsize=(4, 12))
		plt.barh(np.arange(SAMPLE_SIZE), avg_inter_rnap_distance[:SAMPLE_SIZE],
			tick_label=sorted_TU_ids[:SAMPLE_SIZE],
			color=bar_colors)
		plt.xlabel("Average distance between RNAPs (nt)")
		plt.axvline(RNAP_footprint_size, linestyle='--', color='k')

		# Add values to each bar
		for i, v in enumerate(avg_inter_rnap_distance[:SAMPLE_SIZE]):
			plt.text(v - 1, i, "{0:.1f}".format(v), color='white', fontsize=5,
				horizontalalignment='right', verticalalignment='center')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
