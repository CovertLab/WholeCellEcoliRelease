"""
Analysis plot to check the effects of transient gene dosage on transcription
probabilities of RNAs.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/11/18
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


RNA_ID_LIST = ['RRFA-RRNA', 'RRLA-RRNA', 'RRSA-RRNA', # rRNAs
	'alaU-tRNA', # tRNA
	'EG10864_RNA', # ribosomal protein
	'EG10893_RNA', # RNA polymerase
	'EG10375_RNA', 'EG11586_RNA', 'EG11315_RNA',
	'G7559_RNA', 'G7432_RNA', 'RNA0-126', 'G7179_RNA',
	'EG11981_RNA', 'G6985_RNA', 'EG11744_RNA'] # mRNAs, with different positions

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Read from sim_data
		rna_ids = {rna: i for i, rna in enumerate(sim_data.process.transcription.rnaData["id"])}
		rna_idx = [rna_ids[x + "[c]"] for x in RNA_ID_LIST]
		rna_coordinates = sim_data.process.transcription.rnaData[
			"replicationCoordinate"][rna_idx]

		forward_sequence_length = sim_data.process.replication.replichore_lengths[0]
		reverse_sequence_length = sim_data.process.replication.replichore_lengths[1]

		relative_positions = np.array([float(x)/forward_sequence_length
			if x > 0 else float(-x)/reverse_sequence_length
			for x in rna_coordinates])

		all_parca_synth_probs = sim_data.process.transcription.rnaSynthProb[sim_data.condition]

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		synth_prob_reader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))

		# Load data
		time = main_reader.readColumn('time')

		gene_copy_numbers = synth_prob_reader.readColumn(
			"gene_copy_number")[:, rna_idx]
		synth_probs = synth_prob_reader.readColumn(
			"rnaSynthProb")[:, rna_idx]
		parca_synth_probs = all_parca_synth_probs[rna_idx]

		n_plots = len(RNA_ID_LIST)

		fig = plt.figure()
		fig.set_size_inches(8, 3 * n_plots)
		gs = gridspec.GridSpec(n_plots, 1)

		for i, (rna_id, rna_pos) in enumerate(zip(RNA_ID_LIST, relative_positions)):
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
			ax2.set_title("%s, position = %.2f" % (rna_id, rna_pos))
			ax2.plot(time, gene_copy_numbers[:, i], color='r', label="Gene dosage")
			ax2.legend(loc=1)

		fig.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
