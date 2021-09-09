"""
Analysis plot to check the effects of transient gene dosage on transcription
probabilities of RNAs.
"""

from __future__ import absolute_import, division, print_function


import os

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from six.moves import cPickle, zip

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader


RNA_IDS = ['RRFA-RRNA', 'RRLA-RRNA', 'RRSA-RRNA', 'alaU-tRNA']

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Read from sim_data
		rna_replication_coordinates = sim_data.process.transcription.rna_data[
			"replication_coordinate"]
		rna_ids = sim_data.process.transcription.rna_data["id"]
		forward_sequence_length = sim_data.process.replication.replichore_lengths[0]
		reverse_sequence_length = sim_data.process.replication.replichore_lengths[1]

		# Pick representative genes from different locations on the genome
		rna_relative_positions = np.array([float(x)/forward_sequence_length
			if x > 0 else float(-x)/reverse_sequence_length
			for x in rna_replication_coordinates])
		selected_rna_indexes = list(np.argsort(rna_relative_positions)[
			::(len(rna_replication_coordinates)//10)])

		# Add extra genes specified in RNA_IDS
		rna_ids_to_indexes = {
			rna: i for i, rna
			in enumerate(sim_data.process.transcription.rna_data["id"])}
		selected_rna_indexes.extend(
			[rna_ids_to_indexes[rna_id + "[c]"] for rna_id in RNA_IDS]
			)
		selected_rna_indexes = np.array(selected_rna_indexes)

		all_parca_synth_probs = sim_data.process.transcription.rna_synth_prob[sim_data.condition]

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		synth_prob_reader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))

		# Load data
		time = main_reader.readColumn('time')

		gene_copy_numbers = synth_prob_reader.readColumn(
			"gene_copy_number")[:, selected_rna_indexes]
		synth_probs = synth_prob_reader.readColumn(
			"rnaSynthProb")[:, selected_rna_indexes]
		parca_synth_probs = all_parca_synth_probs[selected_rna_indexes]

		n_plots = len(selected_rna_indexes)

		fig = plt.figure()
		fig.set_size_inches(8, 3 * n_plots)
		gs = gridspec.GridSpec(n_plots, 1)

		for (i, rna_index) in enumerate(selected_rna_indexes):
			ax1 = plt.subplot(gs[i, 0])
			ax1.set_ylabel("Transcription probability")
			ax1.plot(time, synth_probs[:, i], label="Transcription probability")
			ax1.axhline(parca_synth_probs[i],
				linestyle="--", color='k', linewidth=3,
				label="Fit transcription probability")
			ax1.legend(loc=2)

			ax2 = ax1.twinx()
			ax2.set_xlabel("Time [s]")
			ax2.set_ylabel("Gene dosage (copy number)")
			ax2.set_ylim([0, 10])
			ax2.set_title(
				"%s, position = %.2f"
				% (rna_ids[rna_index], rna_relative_positions[rna_index]))
			ax2.plot(time, gene_copy_numbers[:, i], color='r', label="Gene dosage")
			ax2.legend(loc=1)

		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
