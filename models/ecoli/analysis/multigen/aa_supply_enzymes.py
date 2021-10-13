"""
Show the concentration of amino acid synthesis enzymes compared to the expected
values from the parca for basal and with amino acid conditions. Can be helpful
in troubleshooting issues with synthesis enzyme production and could be paired
with the aa_supply analysis plot.
"""

import pickle

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import multigenAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import CONC_UNITS
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.utils import units


PLOT_UNITS = units.umol / units.L


def plot(gs, rows, cols, times, with_aa_reference, basal_reference, enzyme_conc,
		 probability_per_amino_acid, offset, label, aa_ids):
	start = times[0, 0]
	end = times[-1, 0]

	for i, (with_aa, basal, sim, prob) in enumerate(
			zip(with_aa_reference, basal_reference, enzyme_conc, probability_per_amino_acid)):
		row = i // cols
		col = i % cols
		ax = plt.subplot(gs[row + offset, col])
		ax_right = plt.twinx(ax)
		mean = sim.mean()

		ax.plot(times, sim, label='Simulation')
		ax.plot([start, end], [mean, mean], 'k--', linewidth=1, label='Simulation, mean')
		ax.plot([start, end], [with_aa, with_aa], 'r--', linewidth=1, label='Expected, with AA')
		ax.plot([start, end], [basal, basal], 'g--', linewidth=1, label='Expected, basal')
		ax_right.plot(times, prob, 'k', alpha=0.5, linewidth=0.5,
					  label='Average RNA synth prob (right)')

		ax.set_title(f'{aa_ids[i]} {label}', fontsize=8)
		if row == rows - 1:
			ax.set_xlabel('Time (min)', fontsize=8)
		if col == 0:
			ax.set_ylabel('Concentration (uM)', fontsize=8)
		if col == cols - 1:
			ax_right.set_ylabel('RNA probability', fontsize=8)

		ax.tick_params(labelsize=6)
		ax_right.tick_params(labelsize=6)

	return ax, ax_right


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		metabolism = sim_data.process.metabolism
		transcription = sim_data.process.transcription
		translation = sim_data.process.translation
		complexation = sim_data.process.complexation

		# Amino acid data
		aa_ids = sim_data.molecule_groups.amino_acids
		n_aas = len(aa_ids)

		# Expected expression from parameter calculations
		enzyme_to_amino_acid_fwd = metabolism.enzyme_to_amino_acid_fwd
		enzyme_to_amino_acid_rev = metabolism.enzyme_to_amino_acid_rev
		with_aa_reference_fwd = metabolism.aa_supply_enzyme_conc_with_aa.asNumber(PLOT_UNITS) @ enzyme_to_amino_acid_fwd
		basal_reference_fwd = metabolism.aa_supply_enzyme_conc_basal.asNumber(PLOT_UNITS) @ enzyme_to_amino_acid_fwd
		with_aa_reference_rev = metabolism.aa_supply_enzyme_conc_with_aa.asNumber(PLOT_UNITS) @ enzyme_to_amino_acid_rev
		basal_reference_rev = metabolism.aa_supply_enzyme_conc_basal.asNumber(PLOT_UNITS) @ enzyme_to_amino_acid_rev

		# Get RNA cistrons associated with each enzyme
		monomer_to_cistron = {m['id']: m['cistron_id'] for m in translation.monomer_data}
		mat_i = []
		mat_j = []
		cistron_mapping_indices = {}
		for j, enzyme in enumerate(metabolism.aa_enzymes):
			for monomer in complexation.get_monomers(enzyme)['subunitIds']:
				cistron = monomer_to_cistron[monomer]
				cistron_mapping_indices[cistron] = cistron_mapping_indices.get(cistron, len(cistron_mapping_indices))
				mat_i.append(cistron_mapping_indices[cistron])
				mat_j.append(j)
		cistron_to_enzyme = np.zeros((np.max(mat_i) + 1, np.max(mat_j) + 1))
		cistron_to_enzyme[mat_i, mat_j] = 1
		cistrons = list(cistron_mapping_indices.keys())
		cistron_synth_prob_data_indices = {cistron: i for i, cistron in enumerate(transcription.cistron_data['id'])}
		enzyme_cistron_indices = np.array([cistron_synth_prob_data_indices[cistron] for cistron in cistrons])

		# Mapping matrix from cistrons to RNAs
		cistron_tu_mapping_matrix = transcription.cistron_tu_mapping_matrix

		# Attenuation information
		attenuated_indices = transcription.attenuated_rna_indices
		n_rnas = len(transcription.rna_data)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
		cell_paths = ap.get_cells()

		# Load data
		times = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True) / 60
		counts_to_mol = (CONC_UNITS * read_stacked_columns(
			cell_paths, 'EnzymeKinetics', 'countsToMolar', remove_first=True)).asNumber(PLOT_UNITS)
		enzyme_counts_fwd = read_stacked_columns(
			cell_paths, 'GrowthLimits', 'aa_supply_enzymes_fwd', remove_first=True)
		enzyme_counts_rev = read_stacked_columns(
			cell_paths, 'GrowthLimits', 'aa_supply_enzymes_rev', remove_first=True)
		probabilities = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'rnaSynthProb', remove_first=True)
		attenuation = read_stacked_columns(
			cell_paths, 'TranscriptElongationListener', 'attenuation_probability', remove_first=True)

		# Calculate derived quantities
		enzyme_conc_fwd = (enzyme_counts_fwd * counts_to_mol).T
		enzyme_conc_rev = (enzyme_counts_rev * counts_to_mol).T
		full_attenuation = np.zeros((attenuation.shape[0], n_rnas))
		full_attenuation[:, attenuated_indices] = attenuation
		no_attenuation_probabilities = 1 - full_attenuation
		probabilities_post_attenuation = probabilities * no_attenuation_probabilities
		rnas_per_amino_acid_fwd = (cistron_to_enzyme @ enzyme_to_amino_acid_fwd).sum(axis=0)
		probability_per_amino_acid_fwd = (
			(probabilities_post_attenuation @ cistron_tu_mapping_matrix.T)[:, enzyme_cistron_indices]
			@ cistron_to_enzyme @ enzyme_to_amino_acid_fwd / rnas_per_amino_acid_fwd).T
		rnas_per_amino_acid_rev = (cistron_to_enzyme @ enzyme_to_amino_acid_rev).sum(axis=0)
		probability_per_amino_acid_rev = (
			(probabilities_post_attenuation @ cistron_tu_mapping_matrix.T)[:, enzyme_cistron_indices]
			@ cistron_to_enzyme @ enzyme_to_amino_acid_rev / rnas_per_amino_acid_rev).T

		# Plot data
		plt.figure(figsize=(16, 24))
		n_subplots = n_aas + 1  # 1 for legend
		rows = int(np.ceil(np.sqrt(n_subplots)))
		cols = int(np.ceil(n_subplots / rows))
		gs = gridspec.GridSpec(nrows=2*rows, ncols=cols)

		## Plot data for each amino acid
		plot(gs, rows, cols, times, with_aa_reference_fwd, basal_reference_fwd,
			 enzyme_conc_fwd, probability_per_amino_acid_fwd, 0, 'forward', aa_ids)
		ax, ax_right = plot(gs, rows, cols, times, with_aa_reference_rev, basal_reference_rev,
			 enzyme_conc_rev, probability_per_amino_acid_rev, rows, 'reverse', aa_ids)

		## Display legend
		handles, labels = ax.get_legend_handles_labels()
		handles_right, labels_right = ax_right.get_legend_handles_labels()
		ax = plt.subplot(gs[-1, -1], frameon=False)
		ax.axis('off')
		ax.legend(handles + handles_right, labels + labels_right, loc='center', frameon=False)

		## Save plot
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
