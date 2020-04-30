"""
Plots heatmap of superhelical densities across the chromosome over time

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/29/20
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath

# For the full heatmap plot, superhelical densities are downsampled by
# averaging across a fixed number of base pairs, to bring down the size of the
# image to reasonable numbers. This variable determines this ratio.
DOWNSAMPLING_RATIO = 100

# Range of relative coordinates for cropped heatmap
CROP_LOWER_BOUND = 240500
CROP_UPPER_BOUND = 246000

# Upper and lower bounds for heatmap (values outside this range are clipped)
HEATMAP_UPPER_BOUND = 0.01
HEATMAP_LOWER_BOUND = -0.01

# Global parameters for all subplots
TICK_ATTRS = {'labelsize': 30, 'length': 10}
LABEL_SIZE = 30


def format_heatmap_subplot(ax, time):
	ax.tick_params(**TICK_ATTRS)
	ax.set_xticks([0, len(time) - 1])
	ax.set_xticklabels([0, '%.1f' % (time[-1] / 60,)])
	ax.set_xlabel('Time (min)', size=LABEL_SIZE)
	ax.set_ylabel('Position on chromosome', size=LABEL_SIZE)

def add_colorbar_subplot(ax, fig, heatmap):
	colorbar = fig.colorbar(heatmap, cax=ax, orientation='horizontal')
	colorbar.set_ticks([HEATMAP_LOWER_BOUND, 0, HEATMAP_UPPER_BOUND])
	colorbar.set_label('Superhelical Density $(L - L_0)/L_0$', size=LABEL_SIZE)
	ax.set_xticklabels(
		['<%g' % (HEATMAP_LOWER_BOUND,), 0, '>%g' % (HEATMAP_UPPER_BOUND,)])
	ax.tick_params(**TICK_ATTRS)


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception('simOutDir does not currently exist as a directory')

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Load parameters
		replichore_lengths = sim_data.process.replication.replichore_lengths
		min_coordinates = -replichore_lengths[1]
		max_coordinates = replichore_lengths[0]
		chromosome_length = max_coordinates - min_coordinates

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		supercoiling_reader = TableReader(os.path.join(simOutDir, 'DnaSupercoiling'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		left_boundary_coordinates = supercoiling_reader.readColumn(
			'segment_left_boundary_coordinates')
		right_boundary_coordinates = supercoiling_reader.readColumn(
			'segment_right_boundary_coordinates')
		superhelical_densities = supercoiling_reader.readColumn(
			'segment_superhelical_densities')

		# Initialize arrays
		downsampled_density_array = np.zeros(
			(chromosome_length // DOWNSAMPLING_RATIO + 1, len(time)))
		cropped_density_array = np.zeros(
			(CROP_UPPER_BOUND - CROP_LOWER_BOUND, len(time)))
		mean_density_time = np.zeros(len(time))
		mean_density_position = np.zeros(chromosome_length)

		# Loop through each timestep
		for i in range(len(time)):
			# Initialize arrays for density and counts of chromosome segments
			# that cover each base pair (to correctly average densities of
			# segments from different chromosomal domains)
			density_this_timestep = np.zeros(chromosome_length, np.float64)
			density_counts = np.zeros(chromosome_length, np.int64)

			# Get data from this timestep
			lb_this_timestep = left_boundary_coordinates[i, :]
			rb_this_timestep = right_boundary_coordinates[i, :]
			sd_this_timestep = superhelical_densities[i, :]

			# Truncate NaNs
			lb_this_timestep = lb_this_timestep[
				np.logical_not(np.isnan(lb_this_timestep))].astype(np.int64)
			rb_this_timestep = rb_this_timestep[
				np.logical_not(np.isnan(rb_this_timestep))].astype(np.int64)
			sd_this_timestep = sd_this_timestep[
				np.logical_not(np.isnan(sd_this_timestep))].astype(np.int64)

			# Loop through each segment
			for lb, rb, sd in zip(lb_this_timestep, rb_this_timestep, sd_this_timestep):
				# Shift coordinates
				index_start = lb - min_coordinates
				index_end = rb - min_coordinates

				# Calculate running average of superhelical densities at each
				# base pair
				density_this_timestep[index_start:index_end] = np.divide(
					np.multiply(
						density_counts[index_start:index_end],
						density_this_timestep[index_start:index_end]
						) + sd,
					density_counts[index_start:index_end] + 1
					)

				# Update counts of segments that cover each base pair
				density_counts[index_start:index_end] += 1

			# Downsample density array by averaging over fixed lengths
			downsampled_density_array[:, i] = np.nanmean(
				np.pad(
					density_this_timestep,
					(0, (DOWNSAMPLING_RATIO - density_this_timestep.size % DOWNSAMPLING_RATIO) % DOWNSAMPLING_RATIO),
					mode='constant', constant_values=np.NaN
					).reshape(-1, DOWNSAMPLING_RATIO),
				axis=1)[::-1]

			# Crop specified region for cropped heatmap
			cropped_density_array[:, i] = density_this_timestep[
				(CROP_LOWER_BOUND - min_coordinates):(CROP_UPPER_BOUND - min_coordinates)][::-1]

			# Calculate the mean density at the current timestep and the
			# running average of each position
			mean_density_time[i] = np.nanmean(density_this_timestep)
			if i == 0:
				mean_density_position = density_this_timestep[::-1]
			else:
				mean_density_position = (
					i * mean_density_position + density_this_timestep[::-1]
					) / (i + 1)

		fig = plt.figure(constrained_layout=True, figsize=(30, 100))
		widths = [1, 10]
		heights = [1, 35, 0.1]
		gs = gridspec.GridSpec(ncols=2, nrows=3,
			width_ratios=widths, height_ratios=heights, figure=fig)

		# Plot main heatmap
		ax_main = fig.add_subplot(gs[1, 1])
		heatmap = ax_main.imshow(
			downsampled_density_array,
			clim=(HEATMAP_LOWER_BOUND, HEATMAP_UPPER_BOUND),
			cmap='coolwarm', aspect='auto')
		ax_main.tick_params(top=True, labeltop=True, right=True, labelright=True)
		ax_main.yaxis.set_label_position('right')
		ax_main.set_yticks(
			[0, replichore_lengths[0] // DOWNSAMPLING_RATIO,
				chromosome_length // DOWNSAMPLING_RATIO])
		ax_main.set_yticklabels(['+terC', 'oriC', '-terC'])
		format_heatmap_subplot(ax_main, time)

		# Add colorbar for main heatmap
		ax_colorbar = fig.add_subplot(gs[2, 1])
		add_colorbar_subplot(ax_colorbar, fig, heatmap)

		# Plot averaged densities of the entire chromosome at each timestep
		ax_time = fig.add_subplot(gs[0, 1])
		ax_time.plot(time / 60., mean_density_time)
		ax_time.set_xlim(0, time[-1] / 60.)
		ax_time.spines['top'].set_visible(False)
		ax_time.spines['right'].set_visible(False)
		ax_time.spines['bottom'].set_visible(False)
		ax_time.tick_params(bottom=False, labelbottom=False, **TICK_ATTRS)
		ax_time.set_ylabel('Sup. Den.', size=LABEL_SIZE)

		# Plot averaged densities across all timesteps for each chromosomal
		# position
		ax_position = fig.add_subplot(gs[1, 0])
		ax_position.plot(mean_density_position, np.arange(chromosome_length)[::-1])
		ax_position.set_ylim(0, chromosome_length - 1)
		ax_position.set_xlim(mean_density_position.max(), mean_density_position.min())  # Flip x-axis
		ax_position.spines['top'].set_visible(False)
		ax_position.spines['right'].set_visible(False)
		ax_position.spines['left'].set_visible(False)
		ax_position.tick_params(left=False, labelleft=False, **TICK_ATTRS)
		ax_position.set_xlabel('Sup. Den.', size=LABEL_SIZE)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# Output cropped version of plot to separate file
		fig_cropped = plt.figure(constrained_layout=True, figsize=(30, 5))
		heights = [12, 1]
		gs = gridspec.GridSpec(
			ncols=1, nrows=2, height_ratios=heights, figure=fig_cropped)

		# Plot cropped heatmap
		ax_cropped = fig_cropped.add_subplot(gs[0, 0])
		heatmap_cropped = ax_cropped.imshow(
			cropped_density_array,
			clim=(HEATMAP_LOWER_BOUND, HEATMAP_UPPER_BOUND),
			cmap='coolwarm', aspect='auto')
		ax_cropped.set_yticks(
			[0, CROP_UPPER_BOUND - CROP_LOWER_BOUND])
		ax_cropped.set_yticklabels(['%d'%(CROP_UPPER_BOUND, ), '%d'%(CROP_LOWER_BOUND, )])
		format_heatmap_subplot(ax_cropped, time)

		# Add colorbar for cropped heatmap
		ax_cropped_colorbar = fig_cropped.add_subplot(gs[1, 0])
		add_colorbar_subplot(ax_cropped_colorbar, fig_cropped, heatmap_cropped)

		exportFigure(plt, plotOutDir, plotOutFileName + '_cropped', metadata)
		plt.close('all')

if __name__ == '__main__':
	Plot().cli()
