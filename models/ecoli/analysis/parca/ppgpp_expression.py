"""
Plots expression expected with ppGpp versus without.

Useful for troublshooting differences in expression with ppGpp.  Low ribosomal
expression can lead to slower growth than expected (points will appear below
diagonal line).
"""

import os
import pickle

import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objs as go

from models.ecoli.analysis import parcaAnalysisPlot


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)

		# Ribosomal genes
		molecule_groups = sim_data.molecule_groups
		protein_30s = molecule_groups.s30_proteins
		protein_50s = molecule_groups.s50_proteins
		rrna = molecule_groups.s30_16s_rRNA + molecule_groups.s50_23s_rRNA + molecule_groups.s50_5s_rRNA
		protein_to_cistron = {p['id']: p['cistron_id'] for p in sim_data.process.translation.monomer_data}
		cistron_ribosome = [protein_to_cistron[p] for p in (protein_30s + protein_50s)]
		rna_ids = sim_data.process.transcription.rna_data['id']
		rna_ribosome = {
			rna_ids[i] for cistron_id in cistron_ribosome
			for i in sim_data.process.transcription.cistron_id_to_rna_indexes(cistron_id)}
		rna_ribosome.update(rrna)
		ribosome_mask = np.array([r in rna_ribosome for r in sim_data.process.transcription.rna_data['id']])

		# Labels for points
		labels = np.array([sim_data.common_names.get_common_name(r[:-3]) for r in sim_data.process.transcription.rna_data['id']])

		# Conditions to calculate expression for - will be individual subplots
		conditions = list(sim_data.condition_active_tfs.keys())
		n_conditions = len(conditions)
		expression = [
			(sim_data.process.transcription.rna_expression[condition], sim_data.calculate_ppgpp_expression(condition))
			for condition in conditions
			]

		# Plotly figure
		n_cols = 2
		titles = [c for c in conditions for _ in range(n_cols)]
		fig = make_subplots(rows=n_conditions, cols=n_cols, subplot_titles=titles)

		## All RNAs
		fig.add_traces([
			go.Scatter(x=x, y=y, mode='markers', text=labels)
			for condition, (x, y) in zip(conditions, expression)],
			rows=list((np.arange(n_conditions)+1).astype(int)), cols=list(np.ones(n_conditions, dtype=int)),
			)
		fig.add_traces([
			go.Scatter(
				x=[min(x.min(), y.min()), max(x.max(), y.max())],
				y=[min(x.min(), y.min()), max(x.max(), y.max())],
				mode='lines')
			for x, y in expression],
			rows=list((np.arange(n_conditions)+1).astype(int)), cols=list(np.ones(n_conditions, dtype=int)),
			)

		## Only rRNA and ribosomal mRNA
		fig.add_traces([
			go.Scatter(x=x[ribosome_mask], y=y[ribosome_mask], mode='markers', text=labels[ribosome_mask])
			for condition, (x, y) in zip(conditions, expression)],
			rows=list((np.arange(n_conditions)+1).astype(int)), cols=list(2*np.ones(n_conditions, dtype=int)),
			)
		fig.add_traces([
			go.Scatter(
				x=[min(x.min(), y.min()), max(x.max(), y.max())],
				y=[min(x.min(), y.min()), max(x.max(), y.max())],
				mode='lines')
			for x, y in expression],
			rows=list((np.arange(n_conditions)+1).astype(int)), cols=list(2*np.ones(n_conditions, dtype=int)),
			)

		## Adjust axes
		for i in range(n_conditions):
			for j in range(n_cols):
				fig.update_xaxes(row=i+1, col=j+1, title='Condition expression')
				fig.update_yaxes(row=i+1, col=j+1, title='ppGpp expression')
		fig.update_layout(showlegend=False)

		# TODO: integrate plotly with exportFigure like matplotlib plt to save pdf and png
		out = os.path.join(plot_out_dir, plot_out_filename)
		fig.write_html(out + '.html')
		# fig.write_html(out + '.pdf')  # Enable pdf if desired


if __name__ == "__main__":
	Plot().cli()
