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
		enzyme_to_amino_acid = metabolism.enzyme_to_amino_acid
		with_aa_reference = metabolism.aa_supply_enzyme_conc_with_aa.asNumber(PLOT_UNITS) @ enzyme_to_amino_acid
		basal_reference = metabolism.aa_supply_enzyme_conc_basal.asNumber(PLOT_UNITS) @ enzyme_to_amino_acid

		# Get RNAs associated with each enzyme
		# TODO (ggsun): This part may need to be adjusted for operons
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
		rna_to_enzyme = np.zeros((np.max(mat_i) + 1, np.max(mat_j) + 1))
		rna_to_enzyme[mat_i, mat_j] = 1
		cistrons = list(cistron_mapping_indices.keys())
		rna_data_indices = {rna: i for i, rna in enumerate(transcription.rna_data['id'])}
		enzyme_rna_indices = np.array([rna_data_indices[cistron + '[c]'] for cistron in cistrons])

		# Attenuation information
		attenuated_indices = transcription.attenuated_rna_indices
		n_rnas = len(transcription.rna_data)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot=True)
		cell_paths = ap.get_cells()

		# Load data
		times = read_stacked_columns(cell_paths, 'Main', 'time', remove_first=True) / 60
		counts_to_mol = (CONC_UNITS * read_stacked_columns(
			cell_paths, 'EnzymeKinetics', 'countsToMolar', remove_first=True)).asNumber(PLOT_UNITS)
		enzyme_counts = read_stacked_columns(
			cell_paths, 'GrowthLimits', 'aa_supply_enzymes', remove_first=True)
		probabilities = read_stacked_columns(
			cell_paths, 'RnaSynthProb', 'rnaSynthProb', remove_first=True)[:, enzyme_rna_indices]
		attenuation = read_stacked_columns(
			cell_paths, 'TranscriptElongationListener', 'attenuation_probability', remove_first=True)

		# Calculate derived quantities
		start = times[0, 0]
		end = times[-1, 0]
		enzyme_conc = (enzyme_counts @ enzyme_to_amino_acid * counts_to_mol).T
		full_attenuation = np.ones((attenuation.shape[0], n_rnas))
		full_attenuation[:, attenuated_indices] = attenuation
		attenuation_probabilities = full_attenuation[:, enzyme_rna_indices]
		rnas_per_amino_acid = (rna_to_enzyme @ enzyme_to_amino_acid).sum(axis=0)
		probability_per_amino_acid = (probabilities * attenuation_probabilities @ rna_to_enzyme @ enzyme_to_amino_acid / rnas_per_amino_acid).T

		# Plot data
		plt.figure(figsize=(16, 12))
		n_subplots = n_aas + 1  # 1 for legend
		rows = int(np.ceil(np.sqrt(n_subplots)))
		cols = int(np.ceil(n_subplots / rows))
		gs = gridspec.GridSpec(nrows=rows, ncols=cols)

		## Plot data for each amino acid
		for i, (with_aa, basal, sim, prob) in enumerate(zip(with_aa_reference, basal_reference, enzyme_conc, probability_per_amino_acid)):
			row = i // cols
			col = i % cols
			ax = plt.subplot(gs[row, col])
			ax_right = plt.twinx(ax)
			mean = sim.mean()

			ax.plot(times, sim, label='Simulation')
			ax.plot([start, end], [mean, mean], 'k--', linewidth=1, label='Simulation, mean')
			ax.plot([start, end], [with_aa, with_aa], 'r--', linewidth=1, label='Expected, with AA')
			ax.plot([start, end], [basal, basal], 'g--', linewidth=1, label='Expected, basal')
			ax_right.plot(times, prob, 'k', alpha=0.5, linewidth=0.5, label='Average RNA synth prob (right)')

			ax.set_title(aa_ids[i], fontsize=8)
			if row == rows - 1:
				ax.set_xlabel('Time (min)', fontsize=8)
			if col == 0:
				ax.set_ylabel('Concentration (uM)', fontsize=8)
			if col == cols - 1:
				ax_right.set_ylabel('RNA probability', fontsize=8)

			ax.tick_params(labelsize=6)
			ax_right.tick_params(labelsize=6)

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
