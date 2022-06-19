"""
Show dynamics of active and inactive transcription factors along with binding
to promoters over time.  Also shows the fractions of active and bound
transcription factors.  Useful for troubleshooting TF dynamics and regulatory
control of transcription in cases where certain proteins might not be produced.
"""

import pickle

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_bulk_molecules, read_stacked_columns


def subplot(gs, legend, title, tf_id, gene_id, t, active, bound, inactive, promoters):
	subplots = gs.subgridspec(nrows=1, ncols=2)

	# Derive values to plot
	total_active = active + bound
	total = total_active if inactive is None else total_active + inactive
	frac_promoters_occupied = bound / promoters
	frac_active = total_active / total
	frac_tf_bound = bound / total_active

	# Plot counts traces
	ax = plt.subplot(subplots[0, 0])
	ax.plot(t, total, alpha=0.5, label='Total TF Counts')
	ax.plot(t, active, alpha=0.5, label='Active/Free TF Counts')
	ax.plot(t, bound, alpha=0.5, label='Bound TF Counts')
	if inactive is not None:
		ax.plot(t, inactive, alpha=0.5, label='Inactive TF Counts')
	ax.set_xlabel('Time (min)', fontsize=10)
	ax.set_ylabel(f'{tf_id} ({gene_id})\nCounts', fontsize=10)
	ax.tick_params(axis='both', labelsize=8)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	handles1, labels1 = ax.get_legend_handles_labels()

	# Plot fraction traces
	ax = plt.subplot(subplots[0, 1])
	ax.plot(t, frac_promoters_occupied, alpha=0.5, label='Fraction promoters bound')
	ax.plot(t, frac_active, alpha=0.5, label='Fraction TFs active')
	ax.plot(t, frac_tf_bound, alpha=0.5, label='Fraction active TFs bound')
	ax.set_xlabel('Time (min)', fontsize=10)
	ax.set_ylabel('Fraction', fontsize=10)
	ax.tick_params(axis='both', labelsize=8)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	handles2, labels2 = ax.get_legend_handles_labels()

	# Handle legends for each column
	legend_gs = legend.subgridspec(nrows=1, ncols=2)
	handles = [handles1, handles2]
	labels = [labels1, labels2]
	for col, (handle, label) in enumerate(zip(handles, labels)):
		ax = plt.subplot(legend_gs[0, col])
		ax.axis('off')
		ax.legend(handle, label, loc='lower center', frameon=False, title=title, title_fontsize=12)


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		t_reg = sim_data.process.transcription_regulation
		tf_ids = t_reg.tf_ids
		tf_type = t_reg.tf_to_tf_type
		active_to_bound = t_reg.active_to_bound
		tf_to_gene_id = t_reg.tf_to_gene_id
		max_type_counts = np.unique(list(tf_type.values()), return_counts=True)[1].max()

		active_ids = []
		inactive_ids = []
		for tf in tf_ids:
			active_ids.append(tf + '[c]')

			if tf_type[tf] == '1CS':
				if tf == active_to_bound[tf]:
					inactive_ids.append(sim_data.process.equilibrium.get_unbound(tf + '[c]'))
				else:
					inactive_ids.append(active_to_bound[tf] + '[c]')
			elif tf_type[tf] == '2CS':
				inactive_ids.append(sim_data.process.two_component_system.active_to_inactive_tf[tf + '[c]'])
			else:
				# Dummy ID for 0CS since there is no inactive form
				inactive_ids.append(tf + '[c]')

		cell_paths = self.ap.get_cells()

		# Read data for all cells
		times = read_stacked_columns(cell_paths, 'Main', 'time') / 60
		bound_counts = read_stacked_columns(cell_paths, 'RnaSynthProb', 'nActualBound')
		promoter_counts = read_stacked_columns(cell_paths, 'RnaSynthProb', 'n_available_promoters')
		active_counts, inactive_counts = read_stacked_bulk_molecules(cell_paths, (active_ids, inactive_ids))

		plt.figure(figsize=(20, 30))

		legend_row = 1
		gs = gridspec.GridSpec(nrows=max_type_counts+legend_row, ncols=3)
		n_0cs = 0
		n_1cs = 0
		n_2cs = 0
		for i, (tf_id, active, bound, inactive, promoters) in enumerate(zip(tf_ids, active_counts.T, bound_counts.T, inactive_counts.T, promoter_counts.T)):
			if tf_type[tf_id] == '0CS':
				inactive = None  # Used dummy counts in arrays above because no inactive form
				n_0cs += 1
				grid = gs[n_0cs, 0]
				legend = gs[0, 0]
			elif tf_type[tf_id] == '1CS':
				n_1cs += 1
				grid = gs[n_1cs, 1]
				legend = gs[0, 1]
			else:
				n_2cs += 1
				grid = gs[n_2cs, 2]
				legend = gs[0, 2]

			subplot(grid, legend, tf_type[tf_id], tf_id, tf_to_gene_id[tf_id],
				times, active, bound, inactive, promoters)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
