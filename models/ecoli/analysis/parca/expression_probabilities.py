"""
Plot to visualize basal probabilities and TF bound delta probabilities for
gene expression.  Useful to see which genes are highly expressed and which
ones have low expression or regulation.  Note that some genes will have 0 basal
probability and some TF-gene pairs will have 0 delta probability, which might
be unexpected.  Also saves tsv files to check which genes and TF-gene pairs
are at 0.
"""

import csv
import os
import pickle

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		t_reg = sim_data.process.transcription_regulation
		transcription = sim_data.process.transcription
		replication = sim_data.process.replication
		tf_to_gene_id = t_reg.tf_to_gene_id
		tf_ids = [tf_to_gene_id[tf] for tf in t_reg.tf_ids]
		basal_prob = t_reg.basal_prob
		reg_i = t_reg.delta_prob['deltaI']
		reg_j = t_reg.delta_prob['deltaJ']
		reg_v = t_reg.delta_prob['deltaV']
		rna_ids = transcription.rna_data['id']
		gene_ids = transcription.rna_data['gene_id']
		gene_to_symbol = {d['name']: d['symbol'] for d in replication.gene_data}
		rna_expression = transcription.rna_expression['basal']

		# Calculate statistics
		n_basal = len(basal_prob)
		n_0_basal = np.sum(basal_prob == 0)
		frac_0_basal = n_0_basal / n_basal
		n_delta = len(reg_v)
		n_0_delta = np.sum(reg_v == 0)
		frac_0_delta = n_0_delta / n_delta

		# Sort regulation by TF
		tf_sort = np.argsort(reg_j)
		tf_sorted_i = reg_i[tf_sort]
		tf_sorted_j = reg_j[tf_sort]
		tf_sorted_v = reg_v[tf_sort]

		plt.figure(figsize=(6, 12))
		gs = gridspec.GridSpec(nrows=3, ncols=1)

		# Plot sorted basal probabilities
		ax = plt.subplot(gs[0, :])
		ax.bar(range(len(basal_prob)), sorted(basal_prob))
		ax.set_yscale('symlog', linthreshy=np.min(basal_prob[basal_prob > 0]))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelsize=6, bottom=False, labelbottom=False)
		ax.set_xlabel('Sorted genes')
		ax.set_ylabel('Basal prob')
		ax.set_title(f'{n_0_basal}/{n_basal} ({100*frac_0_basal:.1f}%) at 0', fontsize=10)

		# Plot delta probabilities grouped by genes
		ax = plt.subplot(gs[1, :])
		ax.bar(range(len(reg_v)), sorted(reg_v))
		ax.set_yscale('symlog', linthreshy=np.min(np.abs(reg_v[np.abs(reg_v) > 0])))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelsize=6, bottom=False, labelbottom=False)
		ax.set_xlabel('Sorted genes')
		ax.set_ylabel('Delta prob')
		ax.set_title(f'{n_0_delta}/{n_delta} ({100*frac_0_delta:.1f}%) at 0', fontsize=10)

		# Plot delta probabilities grouped by transcription factors
		ax = plt.subplot(gs[2, :])
		ticks = np.array([0] + list(np.where(tf_sorted_j[:-1] != tf_sorted_j[1:])[0] + 1) + [len(tf_sorted_j)])
		tf_x = (ticks[1:] + ticks[:-1]) / 2
		for begin, end in zip(ticks[:-1], ticks[1:]):
			ax.bar(range(begin, end), sorted(tf_sorted_v[begin:end]))
		ax.set_yscale('symlog', linthreshy=np.min(np.abs(reg_v[np.abs(reg_v) > 0])))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_xticks(ticks)
		ax.set_xticks(tf_x, minor=True)
		labels = ax.set_xticklabels(tf_ids, minor=True, fontsize=4)
		for i, label in enumerate(labels):
			label.set_y(label.get_position()[1] - (i % 4) * 0.015)
		ax.tick_params(labelsize=6, labelbottom=False)
		ax.tick_params(axis='x', which='minor', bottom=False)
		ax.set_ylabel('Delta prob')

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')

		# Save data to tsv files for easy lookup
		with open(os.path.join(plot_out_dir, f'{plot_out_filename}_basal.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['RNA', 'Gene symbol', 'Basal prob', 'Expression'])
			for idx in np.argsort(basal_prob):
				writer.writerow([
					rna_ids[idx],
					gene_to_symbol[gene_ids[idx]],
					basal_prob[idx],
					rna_expression[idx],
					])

		with open(os.path.join(plot_out_dir, f'{plot_out_filename}_delta.tsv'), 'w') as f:
			writer = csv.writer(f, delimiter='\t')
			writer.writerow(['TF', 'RNA', 'Gene symbol', 'Delta prob'])
			for begin, end in zip(ticks[:-1], ticks[1:]):
				for idx in np.argsort(tf_sorted_v[begin:end]) + begin:
					writer.writerow([
						tf_ids[tf_sorted_j[idx]],
						rna_ids[tf_sorted_i[idx]],
						gene_to_symbol[gene_ids[tf_sorted_i[idx]]],
						tf_sorted_v[idx],
						])

if __name__ == "__main__":
	Plot().cli()
