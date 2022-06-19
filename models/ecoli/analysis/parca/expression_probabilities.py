"""
Plot to visualize basal probabilities and TF bound delta probabilities for
gene expression.  Useful to see which genes are highly expressed and which
ones have low expression or regulation.  Note that some genes will have 0 basal
probability and some TF-gene pairs will have 0 delta probability, which might
be unexpected.  Also saves tsv files to check which genes and TF-gene pairs
are at 0.

TODO:
	include tRNA attenuation basal_prob adjustments and compare?
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

		# sim_data classes used
		t_reg = sim_data.process.transcription_regulation
		transcription = sim_data.process.transcription
		replication = sim_data.process.replication

		## Normal expression
		tf_to_gene_id = t_reg.tf_to_gene_id
		tf_ids = [tf_to_gene_id[tf] for tf in t_reg.tf_ids]
		basal_prob = t_reg.basal_prob
		reg_i = t_reg.delta_prob['deltaI']
		reg_j = t_reg.delta_prob['deltaJ']
		reg_v = t_reg.delta_prob['deltaV']
		rna_ids = transcription.rna_data['id']
		gene_symbols = [
			sim_data.common_names.get_common_name(rna_id[:-3])
			for rna_id in rna_ids]
		rna_expression = transcription.rna_expression['basal']

		## ppGpp expression
		ppgpp_conc = sim_data.growth_rate_parameters.get_ppGpp_conc(sim_data.condition_to_doubling_time['basal'])
		ppgpp_prob, _ = transcription.synth_prob_from_ppgpp(ppgpp_conc, replication.get_average_copy_number)

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

		plt.figure(figsize=(12, 12))
		gs = gridspec.GridSpec(nrows=3, ncols=2)

		# Plot sorted basal probabilities
		min_basal = np.min(basal_prob[basal_prob > 0])
		ax = plt.subplot(gs[0, 0])
		ax.fill_between(range(len(basal_prob)), sorted(basal_prob))
		ax.set_yscale('symlog', linthreshy=min_basal)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelsize=6, bottom=False, labelbottom=False)
		ax.set_xlabel('Sorted genes')
		ax.set_ylabel('Basal prob')
		ax.set_title(f'{n_0_basal}/{n_basal} ({100*frac_0_basal:.1f}%) at 0', fontsize=10)

		# Plot sorted basal probabilities for ppGpp
		min_ppgpp = np.min(ppgpp_prob[ppgpp_prob > 0])
		ax = plt.subplot(gs[1, 0])
		ax.fill_between(range(len(ppgpp_prob)), sorted(ppgpp_prob))
		ax.set_yscale('symlog', linthreshy=min_ppgpp)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelsize=6, bottom=False, labelbottom=False)
		ax.set_xlabel('Sorted genes')
		ax.set_ylabel('ppGpp basal prob')
		ax.set_title(f'{n_0_basal}/{n_basal} ({100*frac_0_basal:.1f}%) at 0', fontsize=10)

		# Compare basal and ppGpp probabilities
		min_val = min(min_basal, min_ppgpp)
		ax = plt.subplot(gs[2, 0])
		ax.plot(basal_prob, ppgpp_prob, 'o', alpha=0.5)
		ax.set_xscale('symlog', linthreshx=min_val)
		ax.set_yscale('symlog', linthreshy=min_val)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_xlabel('Basal prob')
		ax.set_ylabel('ppGpp basal prob')

		# Plot delta probabilities grouped by genes
		ax = plt.subplot(gs[0, 1])
		ax.fill_between(range(len(reg_v)), sorted(reg_v))
		ax.set_yscale('symlog', linthreshy=np.min(np.abs(reg_v[np.abs(reg_v) > 0])))
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.tick_params(labelsize=6, bottom=False, labelbottom=False)
		ax.set_xlabel('Sorted genes')
		ax.set_ylabel('Delta prob')
		ax.set_title(f'{n_0_delta}/{n_delta} ({100*frac_0_delta:.1f}%) at 0', fontsize=10)

		# Plot delta probabilities grouped by transcription factors
		ax = plt.subplot(gs[1, 1])
		ticks = np.array([0] + list(np.where(tf_sorted_j[:-1] != tf_sorted_j[1:])[0] + 1) + [len(tf_sorted_j)])
		tf_x = (ticks[1:] + ticks[:-1]) / 2
		for begin, end in zip(ticks[:-1], ticks[1:]):
			ax.fill_between(range(begin, end), sorted(tf_sorted_v[begin:end]))
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
			writer.writerow(['RNA', 'Gene symbol', 'Basal prob', 'ppGpp basal prob', 'Expression'])
			for idx in np.argsort(basal_prob):
				writer.writerow([
					rna_ids[idx],
					gene_symbols[idx],
					basal_prob[idx],
					ppgpp_prob[idx],
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
						gene_symbols[tf_sorted_i[idx]],
						tf_sorted_v[idx],
						])

if __name__ == "__main__":
	Plot().cli()
