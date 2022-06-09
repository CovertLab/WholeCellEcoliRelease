"""
Template for parca analysis plots
"""

import pickle
import os

from matplotlib import pyplot as plt
# noinspection PyUnresolvedReferences
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from reconstruction.ecoli.dataclasses.process.transcription import RNA_SEQ_ANALYSIS, EXCLUDED_RNA_TYPES
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants


FIGSIZE = (2.5, 3)


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = pickle.load(f)
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Load raw expression levels of individual cistrons from sequencing data
		cistron_raw_expression = []
		cistron_id_to_gene_id = {
			gene['rna_ids'][0]: gene['id'] for gene in raw_data.genes}
		seq_data = {
			x['Gene']: x[sim_data.basal_expression_condition]
			for x in getattr(raw_data.rna_seq_data, f'rnaseq_{RNA_SEQ_ANALYSIS}_mean')}

		gene_id_to_left_end_pos = {
			gene['id']: gene['left_end_pos'] for gene in raw_data.genes
			}
		gene_id_to_right_end_pos = {
			gene['id']: gene['right_end_pos'] for gene in raw_data.genes
			}
		all_cistron_ids = [
			rna['id'] for rna in raw_data.rnas
			if rna['id'] in cistron_id_to_gene_id
			   and gene_id_to_left_end_pos[cistron_id_to_gene_id[rna['id']]] is not None
			   and gene_id_to_right_end_pos[cistron_id_to_gene_id[rna['id']]] is not None
			   and rna['type'] not in EXCLUDED_RNA_TYPES
			]

		for cistron_id in all_cistron_ids:
			gene_id = cistron_id_to_gene_id[cistron_id]
			# If sequencing data is not found, initialize expression to zero.
			cistron_raw_expression.append(seq_data.get(gene_id, 0.))

		cistron_raw_expression = np.array(cistron_raw_expression)

		# Load adjusted expression levels of cistrons from sim_data
		cistron_adjusted_expression = sim_data.process.transcription.cistron_expression['basal']

		# Load mask for adjusted cistrons
		is_adjusted = sim_data.process.transcription.cistron_data['uses_corrected_seq_counts']

		# Rescale adjusted expression to match values with raw expression except
		# for adjusted genes
		index = np.where(~is_adjusted)[0][0]
		adjustment_factor = cistron_raw_expression[index]/cistron_adjusted_expression[index]
		cistron_adjusted_expression *= adjustment_factor

		# Get adjusted values
		before_adjustment = cistron_raw_expression[is_adjusted]
		after_adjustment = cistron_adjusted_expression[is_adjusted]

		fig = plt.figure(figsize=FIGSIZE)

		ax = fig.add_subplot(1, 1, 1)
		ax.scatter(
			np.linspace(-0.25, 0.25, len(before_adjustment)),
			before_adjustment,
			c='#555555', marker='x', linewidth=0.5, s=20, clip_on=False)
		ax.scatter(
			np.ones(len(after_adjustment)),
			after_adjustment,
			c='#555555', marker='x', linewidth=0.5, s=20, clip_on=False)

		ax.set_ylabel('RNAseq counts')
		ax.set_xticks([0, 1])
		ax.set_xticklabels(['before\ncorrection', 'after\ncorrection'])
		ax.set_yticks([0, 500, 1000])
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["bottom"].set_position(("outward", 15))
		ax.spines["left"].set_position(("outward", 30))
		ax.set_xlim([0, 1])
		ax.set_ylim([0, 1000])

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
