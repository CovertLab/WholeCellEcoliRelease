"""
Plots the histograms of the copy number of each RNA at each generation for
multiple-seed simulations.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/2018
"""

from __future__ import absolute_import
from __future__ import division

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import izip, cycle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils.filepath import makedirs
from models.ecoli.analysis import cohortAnalysisPlot

# Number of RNAs sampled for Plot 1
RNA_SAMPLE_COUNT = 10


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Check if the given variant directory exists
		if not os.path.isdir(variantDir):
			raise Exception, "variantDir does not currently exist as a directory."

		# Make plotOut directory if none exists
		makedirs(plotOutDir)

		# Get paths for all cell simulations in each seed
		ap = AnalysisPaths(variantDir, cohort_plot = True)
		n_seed = ap.n_seed
		n_generation = ap.n_generation

		# If the simulation does not have multiple seeds, skip analysis
		if n_seed <= 1:
			print "Skipping -- rnaCopyNumberDistribution only runs for simulations with multiple seeds."
			return

		# Group simulations by generation
		sim_dirs_grouped_by_gen = []
		for gen_idx in range(n_generation):
			sim_dirs_grouped_by_gen.append(ap.get_cells(generation = [gen_idx]))

		# Load simDataFile
		simData = cPickle.load(open(simDataFile))

		# Get IDs for RNA from simData
		ids_rna = simData.process.transcription.rnaData["id"]
		n_rnas = len(ids_rna)

		# Get cell density constant
		cell_density = simData.constants.cellDensity

		# Load simData from first simulation to extract indices
		simOutDir = os.path.join(sim_dirs_grouped_by_gen[0][0], "simOut")
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

		# Get indices for RNA in bulkMolecules
		molecule_ids = bulkMolecules.readAttribute("objectNames")
		idx_rna = [molecule_ids.index(x) for x in ids_rna]

		# Get mass data and calculate initial cell volume (used as standard volume
		# when normalizing RNA counts)
		massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massDataFile.readColumn("cellMass")*units.fg
		expected_initial_volume = cellMass[0]/cell_density

		# Seed np.random
		np.random.seed(21)

		# Extract RNA counts from all simData
		rna_counts = np.zeros((n_generation, n_seed, n_rnas), dtype=np.int)

		# For each generation and seed
		for gen_idx, simDirs in enumerate(sim_dirs_grouped_by_gen):
			for seed_idx, simDir in enumerate(simDirs):

				# Read simData
				simOutDir = os.path.join(simDir, "simOut")

				# Get BulkMolecule and time data
				bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

				# Pick out random timepoint
				idx_timepoint = np.random.randint(0, high=len(time))

				# Get mass data
				massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = massDataFile.readColumn("cellMass")*units.fg

				# Calculate cell volume
				cell_volume = cellMass[idx_timepoint]/cell_density

				# Read counts of all bulk molecules
				bulkCounts = bulkMolecules.readColumn("counts")

				# Read bulkCounts at selected timepoint, and normalize
				bulkCounts_normalized = bulkCounts[idx_timepoint, :]*(expected_initial_volume/cell_volume)

				# Get RNA counts
				rna_counts[gen_idx, seed_idx, :] = bulkCounts_normalized[idx_rna]

		# Calculate statistics
		rna_counts_mean_over_seed = rna_counts.mean(axis=1)
		rna_counts_var_over_seed = rna_counts.var(axis=1)
		rna_counts_std_over_seed = np.sqrt(rna_counts_var_over_seed)
		rna_counts_mean_over_seed_gen = rna_counts_mean_over_seed.mean(axis=0)
		rna_counts_noise_over_seed = rna_counts_var_over_seed/(rna_counts_mean_over_seed**2)
		rna_counts_fano_over_seed = rna_counts_var_over_seed/rna_counts_mean_over_seed

		# Plot 1: Histograms of average RNA counts and individual RNA
		# counts for each generation
		fig = plt.figure()
		fig.set_size_inches(4*n_generation, 4*(RNA_SAMPLE_COUNT + 1))
		gs = gridspec.GridSpec(RNA_SAMPLE_COUNT + 1, n_generation)
		tick_params_plot1 = {"which": "both", "direction": "out", "top": False, "right": False}

		# 1-1: Plot histogram of average RNA counts for each generation
		for gen_idx in range(n_generation):
			ax = plt.subplot(gs[0, gen_idx])
			mean_counts = rna_counts_mean_over_seed[gen_idx, :]
			ax.hist(mean_counts, bins=np.logspace(-3, 6, 19))
			ax.set_title("Generation %d"%(gen_idx,))
			ax.set_xlabel("Mean RNA counts")
			ax.set_xscale("log")
			ax.tick_params(**tick_params_plot1)
			if gen_idx == 0:
				ax.set_ylabel("Number of RNAs")

		# 1-2: Plot histogram of RNA counts for sampled subset of RNAs
		# Sort RNAs based on mean counts, and sample RNAs uniformly
		mean_count_rank = np.argsort(rna_counts_mean_over_seed_gen)[::-1]
		sampled_ranks = np.linspace(0, n_rnas - 1, RNA_SAMPLE_COUNT).astype(int)

		# Plot histogram for each sampled RNA
		for i, rank in enumerate(sampled_ranks):
			idx_rna_sampled = mean_count_rank[rank]

			# Determine histogram range based on maximum count of the specific RNA
			max_count = rna_counts[:, :, idx_rna_sampled].max()
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
				ax = plt.subplot(gs[i + 1, gen_idx])
				seed_counts = rna_counts[gen_idx, :, idx_rna_sampled]

				# The weights rescale histogram such that all columns sum to one
				weights = np.ones_like(seed_counts)/float(len(seed_counts))

				binProb, _, _ = ax.hist(seed_counts, bins=bins, range=hist_range, color='m', weights=weights)
				max_bin_prob_local = binProb.max()

				# Update maximum bin probability
				if max_bin_prob_local > max_bin_prob:
					max_bin_prob = max_bin_prob_local

				ax.set_xlabel(ids_rna[idx_rna_sampled])

				if max_count < 10:
					ax.set_xlim([0, 10])

				ax.tick_params(**tick_params_plot1)

				if gen_idx == 0:
					ax.set_ylabel("Bin probability")

			# Go back and reset y axis upper limit
			for gen_idx in range(n_generation):
				ax = plt.subplot(gs[i + 1, gen_idx])
				ax.set_ylim([0, 1.1*max_bin_prob])

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_histograms", metadata)
		plt.close()

		# Plot 2: Log fold-changes of means and stds of RNA counts with generation
		fig = plt.figure()
		fig.set_size_inches(8, 10)
		gs = gridspec.GridSpec(2, 1)

		# Calculate log fold-change with respect to first generation
		# Note: RNAs with mean and stds of zero have to be filtered out
		# to avoid divide-by-zero
		rna_counts_mean_over_seed_nonzero = rna_counts_mean_over_seed[:, ~np.any(rna_counts_mean_over_seed == 0, axis=0)]
		rna_counts_std_over_seed_nonzero = rna_counts_std_over_seed[:, ~np.any(rna_counts_std_over_seed == 0, axis=0)]
		rna_counts_mean_over_seed_logfoldchange = np.log10(rna_counts_mean_over_seed_nonzero/rna_counts_mean_over_seed_nonzero[0,:])
		rna_counts_std_over_seed_logfoldchange = np.log10(rna_counts_std_over_seed_nonzero/rna_counts_std_over_seed_nonzero[0,:])

		# Select a subset RNAs to plot for readability
		plot_every_n = 20

		ax1 = plt.subplot(gs[0, 0])
		ax1.plot(np.arange(n_generation), rna_counts_mean_over_seed_logfoldchange[:, ::plot_every_n])
		ax1.plot(np.arange(n_generation), np.zeros(n_generation), linestyle='--', linewidth=5, color='k')
		ax1.set_xlabel("Generation #")
		ax1.set_xticks(np.arange(0, n_generation, step=1))
		ax1.set_ylabel("Log10 fold change, mean RNA counts")
		ax1.set_ylim([-1, 1])

		ax2 = plt.subplot(gs[1, 0])
		ax2.plot(np.arange(n_generation), rna_counts_std_over_seed_logfoldchange[:, ::plot_every_n])
		ax2.plot(np.arange(n_generation), np.zeros(n_generation), linestyle='--', linewidth=5, color='k')
		ax2.set_xlabel("Generation #")
		ax2.set_xticks(np.arange(0, n_generation, step=1))
		ax2.set_ylabel("Log10 fold change, std RNA counts")
		ax2.set_ylim([-1, 1])

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_logFoldChange", metadata)
		plt.close()

		# Plot 3: Plot RNA count noise vs. RNA average copy number
		fig = plt.figure()
		fig.set_size_inches(10, 9*(n_generation + 1))
		gs = gridspec.GridSpec(n_generation + 1, 1)
		xlim_plot3 = [1e-3, 1e4]
		ylim_plot3 = [1e-3, 1e4]

		# Plot for each generation
		for (gen_idx, color) in izip(range(n_generation), cycle('bmyk')):
			ax = plt.subplot(gs[gen_idx, 0])
			ax.scatter(rna_counts_mean_over_seed[gen_idx, :], rna_counts_noise_over_seed[gen_idx, :], s=10, color=color, marker='o', lw=0)
			ax.set_xlabel(r"Mean RNA count ($\mu$)")
			ax.set_xscale("log")
			ax.set_xlim(xlim_plot3)
			ax.set_ylabel(r"RNA noise ($\sigma^2/\mu^2$)")
			ax.set_yscale("log")
			ax.set_ylim(ylim_plot3)
			ax.set_title("Generation %d" % (gen_idx,))
			ax.grid(True)

		# Compare plots between first and last generations
		ax = plt.subplot(gs[-1, 0])
		for (gen_idx, color) in izip([0, n_generation - 1], ['b', 'k']):
			ax.scatter(rna_counts_mean_over_seed[gen_idx, :], rna_counts_noise_over_seed[gen_idx, :],
				s=10, alpha=0.3, color=color, marker='o', label="Generation %d" % (gen_idx,), lw=0)

		ax.set_xlabel(r"Mean RNA count ($\mu$)")
		ax.set_xscale("log")
		ax.set_xlim(xlim_plot3)
		ax.set_ylabel(r"RNA noise ($\sigma^2/\mu^2$)")
		ax.set_yscale("log")
		ax.set_ylim(ylim_plot3)
		ax.legend()
		ax.grid(True)

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_RNANoise", metadata)
		plt.close()

		# Plot 4: Plot RNA Fano factor
		fig = plt.figure()
		fig.set_size_inches(10, 18*n_generation)
		gs = gridspec.GridSpec(2*n_generation, 1)
		fano_upperlim = 5

		# 4-1: Plot Fano factor vs. average copy number
		for gen_idx in range(n_generation):
			ax = plt.subplot(gs[gen_idx, 0])
			ax.scatter(rna_counts_mean_over_seed[gen_idx, :], rna_counts_fano_over_seed[gen_idx, :], s=10, marker='o', lw=0)
			ax.set_xlabel("Mean RNA count ($\mu$)")
			ax.set_xscale("log")
			ax.set_xlim([1e-3, 1e4])
			ax.set_ylabel("RNA Fano Factor ($\sigma^2/\mu$)")
			ax.set_ylim([0, fano_upperlim])
			ax.set_title("Generation %d" % (gen_idx,))
			ax.grid(True)

		# 4-2: Plot histogram for Fano factor
		for gen_idx in range(n_generation):
			ax = plt.subplot(gs[gen_idx + n_generation, 0])
			ax.hist(rna_counts_fano_over_seed[gen_idx, :], bins=50, range=(0, fano_upperlim))
			ax.set_xlabel("RNA Fano Factor ($\sigma^2/\mu$)")
			ax.set_title("Generation %d" % (gen_idx,))

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName + "_FanoFactor", metadata)
		plt.close()


if __name__ == "__main__":
	Plot().cli()
