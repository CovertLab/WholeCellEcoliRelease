"""
Plots the histograms of the copy number of each protein at each generation for
multiple-seed simulations.
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import cycle
from six.moves import cPickle, range, zip

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import units

# Number of proteins sampled for Plot 1
PROTEIN_SAMPLE_COUNT = 10


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get paths for all cell simulations in each seed
		n_seed = self.ap.n_seed
		n_generation = self.ap.n_generation

		# If the simulation does not have multiple seeds, skip analysis
		if n_seed <= 1:
			print("Skipping -- proteinCopyNumberDistribution only runs for simulations with multiple seeds.")
			return

		# Divide simulations by generation number
		sim_dirs_grouped_by_gen = []
		for gen_idx in range(n_generation):
			sim_dirs_grouped_by_gen.append(self.ap.get_cells(generation = [gen_idx]))

		# Load simDataFile and get constants
		sim_data = cPickle.load(open(simDataFile, 'rb'))
		cell_density = sim_data.constants.cell_density
		ids_translation = sim_data.process.translation.monomer_data["id"]
		n_proteins = len(ids_translation)

		# Load simData from first simulation
		simOutDir = os.path.join(sim_dirs_grouped_by_gen[0][0], "simOut")

		# Get mass data and calculate initial cell volume (used as standard volume
		# when normalizing protein counts)
		massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massDataFile.readColumn("cellMass")*units.fg
		expected_initial_volume = cellMass[0]/cell_density

		# Seed np.random
		np.random.seed(21)

		# Extract protein counts from all simData
		protein_counts = np.zeros((n_generation, n_seed, n_proteins), dtype=int)

		for gen_idx, simDirs in enumerate(sim_dirs_grouped_by_gen):
			for seed_idx, simDir in enumerate(simDirs):

				simOutDir = os.path.join(simDir, "simOut")

				# Read counts of protein monomers at chosen timepoint
				monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
				counts = monomerCounts.readColumn("monomerCounts")
				time = monomerCounts.readColumn("time")

				# Pick out random timepoint
				idx_timepoint = np.random.randint(0, high=len(time))

				# Get mass data
				massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = massDataFile.readColumn("cellMass")*units.fg

				# Calculate cell volume chosen timepoint
				cell_volume = cellMass[idx_timepoint]/cell_density

				# Read bulkCounts at selected timepoint, and normalize
				normalized_counts = counts[idx_timepoint, :]*(expected_initial_volume/cell_volume)

				# Get protein monomer counts for calculations now that all complexes are dissociated
				protein_counts[gen_idx, seed_idx, :] = normalized_counts

		# Calculate statistics
		protein_counts_mean_over_seed = protein_counts.mean(axis=1)
		protein_counts_var_over_seed = protein_counts.var(axis=1)
		protein_counts_std_over_seed = np.sqrt(protein_counts_var_over_seed)
		protein_counts_mean_over_seed_gen = protein_counts_mean_over_seed.mean(axis=0)
		protein_counts_noise_over_seed = protein_counts_var_over_seed / (protein_counts_mean_over_seed ** 2)
		protein_counts_fano_over_seed = protein_counts_var_over_seed / protein_counts_mean_over_seed

		# Plot 1: Histograms of average protein counts and individual protein
		# counts for each generation
		fig = plt.figure()
		fig.set_size_inches(4*n_generation, 4*(PROTEIN_SAMPLE_COUNT + 1))
		gs = gridspec.GridSpec(PROTEIN_SAMPLE_COUNT + 1, n_generation)
		tick_params_plot1 = {"which": "both", "direction": "out", "top": False, "right": False}

		# 1-1: Plot histogram of average protein counts for each generation
		for gen_idx in range(n_generation):
			ax = plt.subplot(gs[0, gen_idx])
			mean_counts = protein_counts_mean_over_seed[gen_idx, :]
			ax.hist(mean_counts, bins=np.logspace(-3, 6, 19))
			ax.set_title("Generation %d"%(gen_idx,))
			ax.set_xlabel("Mean protein counts")
			ax.set_xscale("log")
			ax.tick_params(**tick_params_plot1)
			if gen_idx == 0:
				ax.set_ylabel("Number of proteins")

		# 1-2: Plot histogram of protein counts for sampled subset of proteins
		# Sort proteins based on mean counts, and sample proteins uniformly
		mean_count_rank = np.argsort(protein_counts_mean_over_seed_gen)[::-1]
		sampled_ranks = np.linspace(0, n_proteins - 1, PROTEIN_SAMPLE_COUNT).astype(int)

		# Plot histogram for each sampled protein
		for i, rank in enumerate(sampled_ranks):
			idx_sampled_protein = mean_count_rank[rank]

			# Determine histogram range based on maximum count of the specific protein
			max_count = protein_counts[:, :, idx_sampled_protein].max()
			if max_count < 10:
				bins = 10
				hist_range = (0, 10)
			elif max_count < 30:
				bins = max_count
				hist_range = (0, max_count)
			else:
				bins = 15
				hist_range = (0, max_count)

			# Initialize maximum bin probability to zero - will be used to set
			# upper limit for the y axis
			max_bin_prob = 0.0

			# Plot histogram for each generation
			for gen_idx in range(n_generation):
				ax = self.subplot(gs[i + 1, gen_idx])
				seed_counts = protein_counts[gen_idx, :, idx_sampled_protein]

				# The weights rescale histogram such that all columns sum to one
				weights = np.ones_like(seed_counts)/float(len(seed_counts))

				binProb, _, _ = ax.hist(seed_counts, bins=bins, range=hist_range, color='m', weights=weights)
				max_bin_prob_local = binProb.max()

				# Update maximum bin probability
				if max_bin_prob_local > max_bin_prob:
					max_bin_prob = max_bin_prob_local

				ax.set_xlabel(ids_translation[idx_sampled_protein])

				if max_count < 10:
					ax.set_xlim([0, 10])

				ax.tick_params(**tick_params_plot1)

				if gen_idx == 0:
					ax.set_ylabel("Bin probability")

			# Go back and reset y axis upper limit
			for gen_idx in range(n_generation):
				ax = self.subplot(gs[i + 1, gen_idx])
				ax.set_ylim([0, 1.1*max_bin_prob])

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_histograms", metadata)
		plt.close()

		# Plot 2: Log fold-changes of means and stds of protein counts with generation
		fig = plt.figure()
		fig.set_size_inches(8, 10)
		gs = gridspec.GridSpec(2, 1)

		# Calculate log fold-change with respect to first generation
		# Note: proteins with mean and stds of zero have to be filtered out
		# to avoid divide-by-zero
		protein_counts_mean_over_seed_nonzero = protein_counts_mean_over_seed[:, ~np.any(protein_counts_mean_over_seed == 0, axis=0)]
		protein_counts_std_over_seed_nonzero = protein_counts_std_over_seed[:, ~np.any(protein_counts_std_over_seed == 0, axis=0)]
		protein_counts_mean_over_seed_logfoldchange = np.log10(protein_counts_mean_over_seed_nonzero/protein_counts_mean_over_seed_nonzero[0, :])
		protein_counts_std_over_seed_logfoldchange = np.log10(protein_counts_std_over_seed_nonzero/protein_counts_std_over_seed_nonzero[0, :])

		# Select a subset RNAs to plot for readability
		plot_every_n = 20

		ax1 = plt.subplot(gs[0, 0])
		ax1.plot(np.arange(n_generation), protein_counts_mean_over_seed_logfoldchange[:, ::plot_every_n])
		ax1.plot(np.arange(n_generation), np.zeros(n_generation), linestyle='--', linewidth=5, color='k')
		ax1.set_xlabel("Generation #")
		ax1.set_xticks(np.arange(0, n_generation, step=1))
		ax1.set_ylabel("Log10 fold change, mean protein counts")
		ax1.set_ylim([-1, 1])

		ax2 = plt.subplot(gs[1, 0])
		ax2.plot(np.arange(n_generation), protein_counts_std_over_seed_logfoldchange[:, ::plot_every_n])
		ax2.plot(np.arange(n_generation), np.zeros(n_generation), linestyle='--', linewidth=5, color='k')
		ax2.set_xlabel("Generation #")
		ax2.set_xticks(np.arange(0, n_generation, step=1))
		ax2.set_ylabel("Log10 fold change, std protein counts")
		ax2.set_ylim([-1, 1])

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_logFoldChange", metadata)
		plt.close()

		# Plot 3: Plot protein count noise vs. protein average copy number
		fig = plt.figure()
		fig.set_size_inches(10, 9*(n_generation + 1))
		gs = gridspec.GridSpec(n_generation + 1, 1)
		xlim_plot3 = [1e-2, 1e7]
		ylim_plot3 = [1e-5, 1e3]

		# Plot for each generation
		for (gen_idx, color) in zip(range(n_generation), cycle('bmyk')):
			ax = plt.subplot(gs[gen_idx, 0])
			ax.scatter(protein_counts_mean_over_seed[gen_idx, :], protein_counts_noise_over_seed[gen_idx, :], s=10, color=color, marker='o', lw=0)
			ax.set_xlabel(r"Mean protein count ($\mu$)")
			ax.set_xscale("log")
			ax.set_xlim(xlim_plot3)
			ax.set_ylabel(r"Protein noise ($\sigma^2/\mu^2$)")
			ax.set_yscale("log")
			ax.set_ylim(ylim_plot3)
			ax.set_title("Generation %d" % (gen_idx,))
			ax.grid(True)

		# Compare plots between first and last generations
		ax = plt.subplot(gs[-1, 0])
		for (gen_idx, color) in zip([0, n_generation - 1], ['b', 'k']):
			ax.scatter(protein_counts_mean_over_seed[gen_idx, :], protein_counts_noise_over_seed[gen_idx, :],
				s=10, alpha=0.3, color=color, marker='o', label="Generation %d" % (gen_idx,), lw=0)

		ax.set_xlabel(r"Mean protein count ($\mu$)")
		ax.set_xscale("log")
		ax.set_xlim(xlim_plot3)
		ax.set_ylabel(r"Protein noise ($\sigma^2/\mu^2$)")
		ax.set_yscale("log")
		ax.set_ylim(ylim_plot3)
		ax.legend()
		ax.grid(True)

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_proteinNoise", metadata)
		plt.close()

		# Plot 4: Plot protein Fano factor
		fig = plt.figure()
		fig.set_size_inches(10, 18*n_generation)
		gs = gridspec.GridSpec(2*n_generation, 1)
		fano_upperlim = 100

		# 4-1: Plot Fano factor vs. average copy number
		for gen_idx in range(n_generation):
			ax = plt.subplot(gs[gen_idx, 0])
			ax.scatter(protein_counts_mean_over_seed[gen_idx, :], protein_counts_fano_over_seed[gen_idx, :], s=10, marker='o', lw=0)
			ax.set_xlabel("Mean protein count ($\mu$)")
			ax.set_xscale("log")
			ax.set_xlim([1e-3, 1e4])
			ax.set_ylabel("Protein Fano Factor ($\sigma^2/\mu$)")
			ax.set_ylim([0, fano_upperlim])
			ax.set_title("Generation %d" % (gen_idx,))
			ax.grid(True)

		# 4-2: Plot histogram for Fano factor
		for gen_idx in range(n_generation):
			ax = plt.subplot(gs[gen_idx + n_generation, 0])
			ax.hist(protein_counts_fano_over_seed[gen_idx, :], bins=50, range=(0, fano_upperlim))
			ax.set_xlabel("Protein Fano Factor ($\sigma^2/\mu$)")
			ax.set_title("Generation %d" % (gen_idx,))

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_FanoFactor", metadata)
		plt.close()

if __name__ == "__main__":
	Plot().cli()
