"""
Plots the histograms of the copy number of each protein at each generation for
multiple-seed simulations.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/20/2018
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

# Number of proteins sampled for Plot 1
PROTEIN_SAMPLE_COUNT = 50


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
			print "Skipping -- proteinCopyNumberDistribution only runs for simulations with multiple seeds."
			return

		# Divide simulations by generation number
		sim_dirs_grouped_by_gen = []
		for gen_idx in range(n_generation):
			sim_dirs_grouped_by_gen.append(ap.get_cells(generation = [gen_idx]))

		# Load simDataFile
		sim_data = cPickle.load(open(simDataFile))

		# Get IDs for complex of proteins and constituent protein monomers,
		# and IDs for only the complexed proteins
		ids_complexation = sim_data.process.complexation.moleculeNames
		ids_complexation_complex = sim_data.process.complexation.ids_complexes

		# Get IDs for complex of proteins and small molecules, and associated
		# protein monomers and small molecules, and IDs for only the complexes
		ids_equilibrium = sim_data.process.equilibrium.moleculeNames
		ids_equilibrium_complex = sim_data.process.equilibrium.ids_complexes

		# Get IDs for all protein monomers
		ids_translation = sim_data.process.translation.monomerData["id"].tolist()
		n_proteins = len(ids_translation)

		# Get data for proteins in ribosomes
		data_50s = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s50_fullComplex)
		data_30s = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.s30_fullComplex)
		ids_ribosome_subunit = data_50s["subunitIds"].tolist() + data_30s["subunitIds"].tolist()

		# Get data for proteins in RNAP
		data_rnap = sim_data.process.complexation.getMonomers(sim_data.moleculeIds.rnapFull)
		ids_rnap_subunit = data_rnap["subunitIds"].tolist()

		# Get all monomer stoichiometric matrices for protein complexes
		# These matrices will be used to dissociate complexes into its constituent
		# monomer proteins
		complex_stoich = sim_data.process.complexation.stoichMatrixMonomers()
		equilibrium_stoich = sim_data.process.equilibrium.stoichMatrixMonomers()
		ribosome_subunit_stoich = np.hstack((data_50s["subunitStoich"], data_30s["subunitStoich"]))
		rnap_subunit_stoich = data_rnap["subunitStoich"]

		# Get cell density constant
		cell_density = sim_data.constants.cellDensity

		# Load simData from first simulation to extract indices
		simOutDir = os.path.join(sim_dirs_grouped_by_gen[0][0], "simOut")
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

		# Read molecule IDs and cache "ID: index" in a dictionary
		molecule_ids = bulkMolecules.readAttribute("objectNames")
		molecule_dict = {mol: i for i, mol in enumerate(molecule_ids)}

		# Get indices from IDs of protein-protein complexes and constituent monomers
		idx_complexation = np.array([molecule_dict[x] for x in ids_complexation])
		idx_complexation_complex = np.array([molecule_dict[x] for x in ids_complexation_complex])

		# Get indices from IDs of protein-small molecule complexes and constituent monomers
		idx_equilibrium = np.array([molecule_dict[x] for x in ids_equilibrium])
		idx_equilibrium_complex = np.array([molecule_dict[x] for x in ids_equilibrium_complex])  # Only complexes

		# Get indices from IDs of protein monomers, and proteins from ribosomes and RNAPs
		idx_translation = np.array([molecule_dict[x] for x in ids_translation])
		idx_ribosome_subunit = np.array([molecule_dict[x] for x in ids_ribosome_subunit])
		idx_rnap_subunit = np.array([molecule_dict[x] for x in ids_rnap_subunit])

		# Get mass data and calculate initial cell volume (used as standard volume
		# when normalizing protein counts)
		massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = massDataFile.readColumn("cellMass")*units.fg
		expected_initial_volume = cellMass[0]/cell_density

		# Seed np.random
		np.random.seed(21)

		# Extract protein counts from all simData
		protein_counts = np.zeros((n_generation, n_seed, n_proteins), dtype=np.int)

		for gen_idx, simDirs in enumerate(sim_dirs_grouped_by_gen):
			for seed_idx, simDir in enumerate(simDirs):

				simOutDir = os.path.join(simDir, "simOut")

				# Get BulkMolecule and time data
				bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
				time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

				# Pick out random timepoint
				idx_timepoint = np.random.randint(0, high=len(time))

				# Get mass data
				massDataFile = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = massDataFile.readColumn("cellMass")*units.fg

				# Calculate cell volume chosen timepoint
				cell_volume = cellMass[idx_timepoint]/cell_density

				# Read counts of all bulk molecules at chosen timepoint
				bulkCounts = bulkMolecules.readColumn("counts")[idx_timepoint, :]
				bulkCounts = bulkCounts[np.newaxis, :]  # Convert to column vector

				# Load unique molecule data for RNAP and ribosomes
				# Note: Inactive ribosomes and RNA polymerases are accounted for
				# in molecule counts in Complexes, not UniqueMoleculeCounts
				uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
				idx_ribosome = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
				idx_rnap = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
				n_active_ribosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[idx_timepoint, idx_ribosome]
				n_active_rnap = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[idx_timepoint, idx_rnap]

				# Dissociate all complexes
				# Stoichiometric matrices are in dtype float64 - needs casting to integers
				complex_monomer_counts = np.dot(np.negative(bulkCounts[:, idx_complexation_complex]), complex_stoich.T)
				equilibrium_monomer_counts = np.dot(np.negative(bulkCounts[:, idx_equilibrium_complex]), equilibrium_stoich.T)
				bulkCounts[:, idx_complexation] += complex_monomer_counts.astype(np.int)
				bulkCounts[:, idx_equilibrium] += equilibrium_monomer_counts.astype(np.int)

				# Add subunits from RNAP and ribosomes
				n_ribosome_subunit = n_active_ribosome*ribosome_subunit_stoich[np.newaxis, :]
				n_rnap_subunit = n_active_rnap*rnap_subunit_stoich[np.newaxis, :]
				bulkCounts[:, idx_ribosome_subunit] += n_ribosome_subunit.astype(np.int)
				bulkCounts[:, idx_rnap_subunit] += n_rnap_subunit.astype(np.int)

				# Read bulkCounts at selected timepoint, and normalize
				bulkCounts_normalized = bulkCounts*(expected_initial_volume/cell_volume)

				# Get protein monomer counts for calculations now that all complexes are dissociated
				protein_counts[gen_idx, seed_idx, :] = bulkCounts_normalized[:, idx_translation]

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
		fig.set_size_inches(5*n_generation, 5*(PROTEIN_SAMPLE_COUNT + 1))
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
				ax = plt.subplot(gs[i + 1, gen_idx])
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
				ax = plt.subplot(gs[i + 1, gen_idx])
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
		for (gen_idx, color) in izip(range(n_generation), cycle('bmyk')):
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
		for (gen_idx, color) in izip([0, n_generation - 1], ['b', 'k']):
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
