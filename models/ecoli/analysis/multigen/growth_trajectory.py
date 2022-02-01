"""
Show the trajectory of growth rates and cell composition measurements.  Useful
for variants with a shift.

TODO:
	- shares a lot of code with the multigen plot with the same name (make into a function?)
	- add cell division markers on the trace
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns


def plot(ax, x, y, xlabel, ylabel, sim_time=None, timeline=None, ma_time=None):
	# Markers for any timeline shifts
	if sim_time is not None and timeline is not None:
		for t, media in timeline:
			times_after_shift = (sim_time > t)[:len(x)]
			marker_x = x[times_after_shift][:1]
			marker_y = y[times_after_shift][:1]
			ax.plot(marker_x, marker_y, 'rx')
			ax.text(marker_x, marker_y, media, fontsize=6)

	# Settings for different plots and start/end markers for main plot
	ax.plot(x[0], y[0], 'og', markersize=3)
	ax.plot(x[-1], y[-1], 'or', markersize=3)

	# Plot trace and time markers
	trace, = ax.plot(x, y, alpha=0.7)
	if ma_time is not None:
		time_hr = np.floor(ma_time / 3600)
		hour_markers = np.where(np.diff(time_hr))[0] + 1
		ax.plot(x[hour_markers], y[hour_markers], 'o', alpha=0.7, markersize=3, color=trace.get_color())

	# Format axes
	ax.set_xlabel(xlabel, fontsize=8)
	ax.set_ylabel(ylabel, fontsize=8)
	ax.tick_params(labelsize=6)


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		if sim_data.external_state.current_timeline_id:
			timeline = sim_data.external_state.saved_timelines[sim_data.external_state.current_timeline_id]
		else:
			timeline = []

		cell_paths = self.ap.get_cells()

		_, axes = plt.subplots(3, 2, figsize=(8, 12))

		# Load data
		growth_function = lambda x: np.diff(x, axis=0) / x[:-1]
		sim_time = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True).squeeze()
		time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec', remove_first=True).squeeze() / 3600  # hr
		growth = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True).squeeze() * 3600  # 1/hr
		protein = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True).squeeze()
		rna = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True).squeeze()
		protein_growth = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', fun=growth_function).squeeze() / time_step
		rna_growth = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', fun=growth_function).squeeze() / time_step
		small_mol_growth = read_stacked_columns(cell_paths, 'Mass', 'smallMoleculeMass', fun=growth_function).squeeze() / time_step
		growth_means = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate', remove_first=True, fun=np.mean).reshape(-1) * 3600
		protein_means = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True, fun=np.mean).reshape(-1)
		rna_means = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True, fun=np.mean).reshape(-1)
		mass_means = read_stacked_columns(cell_paths, 'Mass', 'cellMass', remove_first=True, fun=np.mean).reshape(-1)

		moving_window = min(201, len(growth))
		convolution_array = np.ones(moving_window) / moving_window

		# Process data
		ratio = rna / protein
		ratio_means = rna_means / protein_means
		growth_ma = np.convolve(growth, convolution_array, mode='valid')
		ratio_ma = np.convolve(ratio, convolution_array, mode='valid')
		protein_growth_ma = np.convolve(protein_growth, convolution_array, mode='valid')
		rna_growth_ma = np.convolve(rna_growth, convolution_array, mode='valid')
		small_mol_growth_ma = np.convolve(small_mol_growth, convolution_array, mode='valid')
		ma_time = np.convolve(sim_time, convolution_array, mode='valid')
		offset_time = sim_time[moving_window-1:]  # offset to get first contribution of time point in ma for timeline shift

		plot(axes[0, 0], mass_means, growth_means,
			 xlabel='Average cell cycle mass (fg)', ylabel='Average cell cycle growth rate (1/hr)')
		plot(axes[1, 0], ratio_means, growth_means,
			 xlabel='Average cell cycle RNA/protein', ylabel='Average cell cycle growth rate (1/hr)')
		plot(axes[2, 0], ratio_ma, growth_ma,
			 xlabel='RNA/protein', ylabel='Growth rate (1/hr)',
			 ma_time=ma_time, sim_time=offset_time, timeline=timeline)
		plot(axes[0, 1], ratio_ma, protein_growth_ma,
			 xlabel='RNA/protein', ylabel='Protein growth rate (1/hr)',
			 ma_time=ma_time, sim_time=offset_time, timeline=timeline)
		plot(axes[1, 1], ratio_ma, rna_growth_ma,
			 xlabel='RNA/protein', ylabel='RNA growth rate (1/hr)',
			 ma_time=ma_time, sim_time=offset_time, timeline=timeline)
		plot(axes[2, 1], ratio_ma, small_mol_growth_ma,
			 xlabel='RNA/protein', ylabel='Small molecule growth rate (1/hr)',
			 ma_time=ma_time, sim_time=offset_time, timeline=timeline)

		for ax in axes.reshape(-1):
			self.remove_border(ax)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
