"""
Generates a table of mRNA transcripts that are present in the model and the
parameters associated with each transcript.
"""

import csv
import pickle
# noinspection PyUnresolvedReferences
import numpy as np
import os

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.utils import units


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		rna_data = sim_data.process.transcription.rna_data
		is_mRNA = rna_data['is_mRNA']
		mRNA_indexes = np.where(is_mRNA)[0]
		mRNA_data = rna_data[is_mRNA]

		mRNA_ids = mRNA_data['id']
		mRNA_lengths = mRNA_data['length'].asNumber(units.nt)
		mRNA_deg_rates = mRNA_data['deg_rate'].asNumber(1/units.s)

		# Normalize expression to sum to 1M
		mRNA_expression = sim_data.process.transcription.rna_expression['basal'][is_mRNA]
		mRNA_expression_norm = mRNA_expression/mRNA_expression.sum() * 10**6

		gene_id_to_name = {
			gene['name']: gene['symbol']
			for gene in sim_data.process.replication.gene_data
			}
		gene_names = [
			gene_id_to_name[x] for x in
			sim_data.process.transcription.cistron_data['gene_id']
		]
		constituent_gene_names = []
		for rna_id in mRNA_ids:
			cistron_indexes = sim_data.process.transcription.rna_id_to_cistron_indexes(rna_id)
			constituent_gene_names.append(
				[gene_names[i] for i in cistron_indexes]
				)

		operon_is_polycistronic = []
		for i in mRNA_indexes:
			for (j, k) in sim_data.process.transcription.operons:
				if i in k:
					operon_is_polycistronic.append(len(j) > 1)
					break

		with open(os.path.join(plot_out_dir, plot_out_filename + '.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow([
				'transcript_id', 'constituent_gene_names', 'length',
				'operon_is_polycistronic', 'normalized_read_counts', 'deg_rate'
				])

			for rna_id, gene_names, length, is_poly, read_counts, deg_rate in zip(
					mRNA_ids, constituent_gene_names, mRNA_lengths,
					operon_is_polycistronic, mRNA_expression_norm,
					mRNA_deg_rates,
					):
				writer.writerow([
					rna_id[:-3],
					', '.join(gene_names),
					length, is_poly, read_counts, deg_rate,
					])


if __name__ == "__main__":
	Plot().cli()
