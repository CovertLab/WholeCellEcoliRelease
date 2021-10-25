"""
Show mean traces and confidence interval for growth related properties over
multiple initial seeds.
"""

import os
import pickle

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def plot_time_series(self, gs, t_flat, y_flat, ylabel, timeline, log_scale=False):
		ax = plt.subplot(gs)

		# Extract y data for each time point (assumes time step lines up across samples)
		data = {}
		for t, y in zip(t_flat, y_flat):
			data.setdefault(t, []).append(y)

		# Calculate mean and standard deviation for each time point
		mean = []
		std = []
		t = np.array(sorted(data.keys()))
		for _t in t:
			d = np.vstack(data[_t])
			mean.append(d.mean(0))
			std.append(d.std(0))
		mean = np.array(mean).squeeze()
		std = np.array(std).squeeze()

		# # Plot all single traces - might be too much for large sets of seeds so commented for now
		# new_cell = np.where(t_flat[:-1] > t_flat[1:])[0] + 1
		# splits = [0] + list(new_cell) + [None]
		# for start, end in zip(splits[:-1], splits[1:]):
		# 	ax.plot(t_flat[start:end], y_flat[start:end], 'k', alpha=0.05, linewidth=0.5)

		# Plot mean as a trace and standard deviation as a shaded area
		ax.plot(t, mean)
		if len(mean.shape) > 1:
			for m, s in zip(mean.T, std.T):
				ax.fill_between(t, m - s, m + s, alpha=0.1)
		else:
			ax.fill_between(t, mean - std, mean + std, alpha=0.1)

		# Format axes
		if log_scale:
			ax.set_yscale('log')
		ax.set_xlabel('Time (min)', fontsize=8)
		ax.set_ylabel(ylabel, fontsize=8)
		ax.tick_params(labelsize=6)
		self.remove_border(ax)

		# Show any media changes
		t_max = t.max()
		y_min, y_max = ax.get_ylim()
		for i, (t_media, media) in enumerate(timeline):
			t_media /= 60
			if t_media < t_max:
				ax.axvline(t_media, color='k', linestyle='--', linewidth=0.5, alpha=0.3)
				if log_scale:
					y_pos = 10**(np.log10(y_max) - (np.log10(y_max) - np.log10(y_min)) * i * 0.03)
				else:
					y_pos = y_max - (y_max - y_min) * i * 0.03
				ax.text(t_media, y_pos, media, fontsize=6)

	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		transcription = sim_data.process.transcription
		ppgpp_id = sim_data.molecule_ids.ppGpp
		rnap_id = sim_data.molecule_ids.full_RNAP
		aa_ids = sim_data.molecule_groups.amino_acids
		ribosome_subunit_ids = [sim_data.molecule_ids.s30_full_complex, sim_data.molecule_ids.s50_full_complex]
		uncharged_trna_names = transcription.rna_data['id'][transcription.rna_data['is_tRNA']]
		charged_trna_names = transcription.charged_trna_names
		aa_from_trna = transcription.aa_from_trna.T
		if sim_data.external_state.current_timeline_id:
			timeline = sim_data.external_state.saved_timelines[
				sim_data.external_state.current_timeline_id]
		else:
			timeline = []

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells()

		# Load attributes
		unique_molecule_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
		unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
		rnap_idx = unique_molecule_ids.index('active_RNAP')
		ribosome_idx = unique_molecule_ids.index('active_ribosome')

		# Load data
		growth_function = lambda x: np.diff(x, axis=0) / x[:-1]
		time = read_stacked_columns(cell_paths, 'Main', 'time',
			remove_first=True, ignore_exception=True).squeeze() / 60
		time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec',
			remove_first=True, ignore_exception=True).squeeze()
		growth_rate = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate',
			remove_first=True, ignore_exception=True).squeeze() * 3600
		protein_growth = read_stacked_columns(cell_paths, 'Mass', 'proteinMass',
			fun=growth_function, ignore_exception=True).squeeze() / time_step * 3600
		rna_growth = read_stacked_columns(cell_paths, 'Mass', 'rnaMass',
			fun=growth_function, ignore_exception=True).squeeze() / time_step * 3600
		small_mol_growth = read_stacked_columns(cell_paths, 'Mass', 'smallMoleculeMass',
			fun=growth_function, ignore_exception=True).squeeze() / time_step * 3600
		ribosome_elong_rate = read_stacked_columns(cell_paths, 'RibosomeData', 'effectiveElongationRate',
			remove_first=True, ignore_exception=True).squeeze()
		rnap_elongations = read_stacked_columns(cell_paths, 'RnapData', 'actualElongations',
			remove_first=True, ignore_exception=True).squeeze()
		counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
			remove_first=True, ignore_exception=True)
		unique_mol_counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
			remove_first=True, ignore_exception=True)
		(ppgpp_counts, uncharged_trna_counts, charged_trna_counts, aa_counts,
			inactive_rnap_counts, ribosome_subunit_counts) = read_stacked_bulk_molecules(cell_paths,
				([ppgpp_id], uncharged_trna_names, charged_trna_names, aa_ids, [rnap_id], ribosome_subunit_ids),
				remove_first=True, ignore_exception=True)

		# Derived quantities
		ppgpp_conc = ppgpp_counts * counts_to_molar.squeeze() * 1000
		charged_trna_counts = charged_trna_counts @ aa_from_trna
		uncharged_trna_counts = uncharged_trna_counts @ aa_from_trna
		fraction_charged = charged_trna_counts / (uncharged_trna_counts + charged_trna_counts)
		aa_conc = counts_to_molar * aa_counts
		active_rnap_counts = unique_mol_counts[:, rnap_idx]
		active_ribosome_counts = unique_mol_counts[:, ribosome_idx]
		rnap_elong_rate = rnap_elongations / time_step / active_rnap_counts
		rnap_fraction_active = active_rnap_counts / (active_rnap_counts + inactive_rnap_counts)
		ribosome_fraction_active = active_ribosome_counts / (active_ribosome_counts + ribosome_subunit_counts.min(1))

		plt.figure(figsize=(15, 15))
		gs = gridspec.GridSpec(4, 3)

		self.plot_time_series(gs[0, 0], time, growth_rate, 'Growth rate\n(1/hr)', timeline)
		self.plot_time_series(gs[1, 0], time, rna_growth, 'RNA growth rate\n(1/hr)', timeline)
		self.plot_time_series(gs[2, 0], time, protein_growth, 'Protein growth rate\n(1/hr)', timeline)
		self.plot_time_series(gs[3, 0], time, small_mol_growth, 'Small mol growth rate\n(1/hr)', timeline)
		self.plot_time_series(gs[0, 1], time, rnap_elong_rate, 'RNAP elongation rate\n(nt/s)', timeline)
		self.plot_time_series(gs[1, 1], time, rnap_fraction_active, 'RNAP active fraction', timeline)
		self.plot_time_series(gs[2, 1], time, ribosome_elong_rate, 'Ribosome elongation rate\n(AA/s)', timeline)
		self.plot_time_series(gs[3, 1], time, ribosome_fraction_active, 'Ribosome active fraction', timeline)
		self.plot_time_series(gs[0, 2], time, ppgpp_conc, 'ppGpp concentration\n(uM)', timeline)
		self.plot_time_series(gs[1, 2], time, fraction_charged, 'Fraction charged', timeline)
		self.plot_time_series(gs[2, 2], time, aa_conc, 'Amino acid concentrations\n(mM)', timeline, log_scale=True)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
