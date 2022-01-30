"""
Show mean traces and confidence interval for growth related properties over
multiple initial seeds.
"""

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.analysis.single.ribosome_limitation import calculate_ribosome_excesses
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


PLOT_SINGLE = False


def trim_axes(axes):
	for ax in axes.flatten():
		for text in ax.texts:
			text.set_visible(False)
	axes[0, 0].set_ylim(0, 2)
	axes[1, 0].set_ylim(0, 2)
	axes[2, 0].set_ylim(0, 2)
	axes[3, 0].set_ylim(0, 2)
	axes[0, 1].set_ylim(0, 100)
	axes[1, 1].set_ylim(0, 1)
	axes[2, 1].set_ylim(0, 25)
	axes[3, 1].set_ylim(0, 1)
	axes[0, 2].set_ylim(0, 1)
	axes[1, 2].set_ylim(0, 1)
	axes[2, 2].set_ylim(0, 10)
	axes[3, 2].set_ylim(0, 1)
	axes[0, 3].set_ylim(0, 300)
	axes[1, 3].set_ylim(0, 1.2)
	axes[2, 3].set_ylim(1e-4, 100)
	axes[3, 3].set_ylim(0, 1)
	axes[0, 4].set_ylim(0, 1)
	axes[1, 4].set_ylim(0, 0.15)
	axes[2, 4].set_ylim(0, 0.3)
	axes[0, 5].set_ylim(0, 10)
	axes[1, 5].set_ylim(0, 0.2)
	axes[2, 5].set_ylim(0, 40)
	axes[3, 5].set_ylim(0, 1)
	axes[0, 6].set_ylim(0, 0.1)
	axes[1, 6].set_ylim(0, 1)
	axes[2, 6].set_ylim(0, 1)
	axes[3, 6].set_ylim(0, 1)
	axes[0, 7].set_ylim(0, 1)
	axes[1, 7].set_ylim(0, 1)
	axes[2, 7].set_ylim(0, 15)
	axes[3, 7].set_ylim(0, 1)


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def plot_time_series(self, ax, t_flat, y_flat, ylabel, timeline, filtered_t, downsample=5, log_scale=False):
		# TODO: add trace labels as arg
		# Extract y data for each time point (assumes time step lines up across samples)
		data = {}
		for t, y in zip(t_flat, y_flat):
			if t not in filtered_t:
				data.setdefault(t, []).append(y)

		# Calculate mean and standard deviation for each group of time points
		# based on downsampling times
		t = []
		mean = []
		std = []
		all_times = np.array(sorted(data.keys()))
		n_dropped = len(all_times) % downsample
		drop_slice = -n_dropped if n_dropped else None
		times = all_times[:drop_slice].reshape(-1, downsample)
		if len(y_flat.shape) > 1:
			stack = np.vstack
		else:
			stack = np.hstack
		for ts in times:
			d = stack([data[_t] for _t in ts])
			t.append(ts.mean())
			mean.append(d.mean(0))
			std.append(d.std(0))
		t = np.array(t)
		mean = np.array(mean).squeeze()
		std = np.array(std).squeeze()

		# Plot all single traces
		if PLOT_SINGLE:
			new_cell = np.where(t_flat[:-1] > t_flat[1:])[0] + 1
			splits = [0] + list(new_cell) + [None]
			for start, end in zip(splits[:-1], splits[1:]):
				ax.plot(t_flat[start:end], y_flat[start:end], 'k', alpha=0.05, linewidth=0.5)

		# Plot mean as a trace and standard deviation as a shaded area
		ax.plot(t, mean)
		if len(mean.shape) > 1:
			for m, s in zip(mean.T, std.T):
				ax.fill_between(t, m - s, m + s, alpha=0.1)
		else:
			ax.axhline(mean.mean(), linestyle='--', color='k', linewidth=0.5)
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

	def plot_hist(self, ax, data, min_val, max_val, label, n_bins=40, sf=2):
		def plot(d):
			patch = ax.hist(d, bins=n_bins, range=(min_val, max_val), alpha=0.7, histtype='step')[-1][0]
			color = patch.get_edgecolor()

			# Add mean +/- std text
			mean = d.mean()
			ax.text(mean, ax.get_ylim()[1], f'{mean:.{sf}f} +/- {d.std():.{sf+1}f}', fontsize=6, color=color)

		if len(data.shape) > 1:
			for d in data.T:
				plot(d)
		else:
			plot(data)

		self.remove_border(ax)
		ax.tick_params(labelsize=6)
		ax.set_xlabel(label, fontsize=8)
		ax.set_ylabel('Count', fontsize=8)


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
		cistron_data = transcription.cistron_data
		if sim_data.external_state.current_timeline_id:
			timeline = sim_data.external_state.saved_timelines[
				sim_data.external_state.current_timeline_id]
		else:
			timeline = []
		rna_fractions = ['is_mRNA', 'is_rRNA', 'is_tRNA']
		convert_to_fraction = lambda x: np.vstack([
			x[:, cistron_data[fraction]].sum(1)
			for fraction in rna_fractions
			]).T
		aa_mw = np.array([sim_data.getter.get_mass(aa[:-3]).asNumber(units.fg / units.count) for aa in aa_ids])
		rna_mw = transcription.rna_data['mw'].asNumber(units.fg / units.count)
		is_mrna = transcription.rna_data['is_mRNA']
		is_rrna = transcription.rna_data['is_rRNA']
		is_trna = transcription.rna_data['is_tRNA']

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells(only_successful=True)

		# Load attributes
		unique_molecule_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
		unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
		rnap_idx = unique_molecule_ids.index('active_RNAP')
		ribosome_idx = unique_molecule_ids.index('active_ribosome')

		# Load data
		growth_function = lambda x: np.diff(x, axis=0) / x[:-1]
		def reduce_rna(mask):
			return lambda x: (x[:, mask] @ rna_mw[mask]).reshape(-1, 1)
		axis_sum = lambda x: x.sum(1).reshape(-1, 1)
		time = read_stacked_columns(cell_paths, 'Main', 'time',
			remove_first=True).squeeze() / 60
		time_step = read_stacked_columns(cell_paths, 'Main', 'timeStepSec',
			remove_first=True).squeeze()
		growth_rate = read_stacked_columns(cell_paths, 'Mass', 'instantaneous_growth_rate',
			remove_first=True).squeeze() * 3600
		protein_growth = read_stacked_columns(cell_paths, 'Mass', 'proteinMass',
			fun=growth_function).squeeze() / time_step * 3600
		rna_growth = read_stacked_columns(cell_paths, 'Mass', 'rnaMass',
			fun=growth_function).squeeze() / time_step * 3600
		small_mol_growth = read_stacked_columns(cell_paths, 'Mass', 'smallMoleculeMass',
			fun=growth_function).squeeze() / time_step * 3600
		protein_mass = read_stacked_columns(cell_paths, 'Mass', 'proteinMass', remove_first=True).squeeze()
		rna_mass = read_stacked_columns(cell_paths, 'Mass', 'rnaMass', remove_first=True).squeeze()
		mrna_mass = read_stacked_columns(cell_paths, 'Mass', 'mRnaMass', remove_first=True).squeeze()
		rrna_mass = read_stacked_columns(cell_paths, 'Mass', 'rRnaMass', remove_first=True).squeeze()
		trna_mass = read_stacked_columns(cell_paths, 'Mass', 'tRnaMass', remove_first=True).squeeze()
		cell_mass = read_stacked_columns(cell_paths, 'Mass', 'cellMass', remove_first=True).squeeze()
		rna_deg_mass = read_stacked_columns(cell_paths, 'RnaDegradationListener', 'countRnaDegraded',
			remove_first=True, fun=reduce_rna(slice(None))).squeeze()
		mrna_deg_mass = read_stacked_columns(cell_paths, 'RnaDegradationListener', 'countRnaDegraded',
			remove_first=True, fun=reduce_rna(is_mrna)).squeeze()
		rrna_deg_mass = read_stacked_columns(cell_paths, 'RnaDegradationListener', 'countRnaDegraded',
			remove_first=True, fun=reduce_rna(is_rrna)).squeeze()
		trna_deg_mass = read_stacked_columns(cell_paths, 'RnaDegradationListener', 'countRnaDegraded',
			remove_first=True, fun=reduce_rna(is_trna)).squeeze()
		ribosome_elong_rate = read_stacked_columns(cell_paths, 'RibosomeData', 'effectiveElongationRate',
			remove_first=True).squeeze()
		rnap_elongations = read_stacked_columns(cell_paths, 'RnapData', 'actualElongations',
			remove_first=True).squeeze()
		aas_elongated = read_stacked_columns(cell_paths, 'GrowthLimits', 'aasUsed', remove_first=True)
		ntps_elongated = read_stacked_columns(cell_paths, 'GrowthLimits', 'ntpUsed', remove_first=True)
		counts_to_molar = read_stacked_columns(cell_paths, 'EnzymeKinetics', 'countsToMolar',
			remove_first=True)
		unique_mol_counts = read_stacked_columns(cell_paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
			remove_first=True)
		rna_fraction_prob = read_stacked_columns(cell_paths, 'RnaSynthProb', 'rna_synth_prob_per_cistron',
			remove_first=True, fun=convert_to_fraction)
		(ppgpp_counts, uncharged_trna_counts, charged_trna_counts, aa_counts,
			inactive_rnap_counts, ribosome_subunit_counts) = read_stacked_bulk_molecules(cell_paths,
				([ppgpp_id], uncharged_trna_names, charged_trna_names, aa_ids, [rnap_id], ribosome_subunit_ids),
				remove_first=True)
		excess, synth_fractions, _ = calculate_ribosome_excesses(sim_data, cell_paths)
		mrna_counts = read_stacked_columns(cell_paths, 'mRNACounts', 'mRNA_counts', remove_first=True, fun=axis_sum).squeeze()
		rna_produced_mass = read_stacked_columns(cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized',
			remove_first=True, fun=reduce_rna(slice(None)))
		mrna_produced_mass = read_stacked_columns(cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized',
			remove_first=True, fun=reduce_rna(is_mrna))
		rrna_produced_mass = read_stacked_columns(cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized',
			remove_first=True, fun=reduce_rna(is_rrna))
		trna_produced_mass = read_stacked_columns(cell_paths, 'TranscriptElongationListener', 'countRnaSynthesized',
			remove_first=True, fun=reduce_rna(is_trna))

		# Derived quantities
		counts_to_molar_squeezed = counts_to_molar.squeeze()
		ppgpp_conc = ppgpp_counts * counts_to_molar_squeezed * 1000
		charged_trna_counts = charged_trna_counts @ aa_from_trna
		uncharged_trna_counts = uncharged_trna_counts @ aa_from_trna
		fraction_charged = charged_trna_counts / (uncharged_trna_counts + charged_trna_counts)
		aa_conc = counts_to_molar * aa_counts
		active_rnap_counts = unique_mol_counts[:, rnap_idx]
		active_ribosome_counts = unique_mol_counts[:, ribosome_idx]
		rnap_elong_rate = rnap_elongations / time_step / active_rnap_counts
		rnap_fraction_active = active_rnap_counts / (active_rnap_counts + inactive_rnap_counts)
		inactive_ribosome_counts = ribosome_subunit_counts.min(1)
		ribosome_fraction_active = active_ribosome_counts / (active_ribosome_counts + inactive_ribosome_counts)
		rnap_conc = counts_to_molar_squeezed * (active_rnap_counts + inactive_rnap_counts) * 1000
		ribosome_conc = counts_to_molar_squeezed * (active_ribosome_counts + inactive_ribosome_counts) * 1000
		rnap_output = counts_to_molar_squeezed * ntps_elongated.sum(1) / time_step
		ribosome_output = counts_to_molar_squeezed * aas_elongated.sum(1) / time_step
		rp_ratio = rna_mass / protein_mass
		aa_mass = aa_counts @ aa_mw
		rpa_ratio = rna_mass / (protein_mass + aa_mass)
		ribosome_fraction_active = active_ribosome_counts / (active_ribosome_counts + ribosome_subunit_counts.min(1))
		protein_fraction = protein_mass / cell_mass
		rna_fraction = rna_mass / cell_mass
		aa_fraction = (protein_mass + aa_mass) / cell_mass
		mrna_rrna_ratio = mrna_mass / rrna_mass
		mrna_fraction = mrna_mass / rna_mass
		rrna_fraction = rrna_mass / rna_mass
		trna_fraction = trna_mass / rna_mass
		mrna_conc = mrna_counts / cell_mass
		rna_deg_rate = rna_deg_mass / rna_mass / time_step * 3600
		mrna_deg_ratio = mrna_deg_mass / rna_deg_mass
		rrna_deg_ratio = rrna_deg_mass / rna_deg_mass
		trna_deg_ratio = trna_deg_mass / rna_deg_mass
		mrna_fraction_produced = mrna_produced_mass / rna_produced_mass
		rrna_fraction_produced = rrna_produced_mass / rna_produced_mass
		trna_fraction_produced = trna_produced_mass / rna_produced_mass
		unique_time, cell_count = np.unique(time, return_counts=True)
		filtered = set(unique_time[cell_count != cell_count.max()])

		def subplots(filtered, downsample=5):
			_, axes = plt.subplots(4, 8, figsize=(25, 15))
			self.plot_time_series(axes[0, 0], time, growth_rate, 'Growth rate\n(1/hr)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 0], time, rna_growth, 'RNA growth rate\n(1/hr)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 0], time, protein_growth, 'Protein growth rate\n(1/hr)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 0], time, small_mol_growth, 'Small mol growth rate\n(1/hr)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[0, 1], time, rnap_elong_rate, 'RNAP elongation rate\n(nt/s)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 1], time, rnap_fraction_active, 'RNAP active fraction',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 1], time, ribosome_elong_rate, 'Ribosome elongation rate\n(AA/s)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 1], time, ribosome_fraction_active, 'Ribosome active fraction',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[0, 2], time, fraction_charged[:, 0], f'Fraction charged\n{aa_ids[0][:-3]} tRNA',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 2], time, fraction_charged[:, 10], f'Fraction charged\n{aa_ids[10][:-3]} tRNA',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 2], time, aa_conc[:, 0], f'{aa_ids[0][:-3]} concentration\n(mM)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 2], time, aa_conc[:, 10], f'{aa_ids[10][:-3]} concentration\n(mM)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[0, 3], time, ppgpp_conc, 'ppGpp concentration\n(uM)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 3], time, fraction_charged, 'Fraction charged',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 3], time, aa_conc, 'Amino acid concentrations\n(mM)',
				timeline, filtered, downsample=downsample, log_scale=True)
			self.plot_time_series(axes[3, 3], time, rna_fraction_prob, 'RNA fraction\nsynthesis probability',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[0, 4], time, rp_ratio, '',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[0, 4], time, rpa_ratio, 'RNA/protein mass fraction\n(with and without free AA)',
				[], filtered, downsample=downsample)
			self.plot_time_series(axes[1, 4], time, rna_fraction, 'RNA mass fraction',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 4], time, protein_fraction, '',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 4], time, aa_fraction, 'Protein mass fraction\n(with and without free AA)',
				[], filtered, downsample=downsample)
			self.plot_time_series(axes[3, 4], unique_time, cell_count, '# cells',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[0, 5], time, rnap_conc, 'RNAP conc\n(uM)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 5], time, rnap_output, 'RNAP output\n(mM NTPs/s)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 5], time, ribosome_conc, 'Ribosome conc\n(uM)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 5], time, ribosome_output, 'Ribosome output\n(mM AA/s)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[0, 6], time, mrna_rrna_ratio, 'mRNA:rRNA ratio',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 6], time, mrna_fraction, '',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 6], time, rrna_fraction, '',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 6], time, trna_fraction, 'RNA mass fractions',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 6], time, rna_deg_rate, 'RNA deg rate',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 6], time, mrna_deg_ratio, '',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 6], time, rrna_deg_ratio, '',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 6], time, trna_deg_ratio, 'RNA deg ratio',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[0, 7], time, excess, 'Excess ribosome RNA/protein',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[1, 7], time, synth_fractions, 'Synthesis fraction RNA/protein/enzymes',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[2, 7], time, mrna_conc, 'mRNA conc\n(count/fg)',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 7], time, mrna_fraction_produced, '',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 7], time, rrna_fraction_produced, '',
				timeline, filtered, downsample=downsample)
			self.plot_time_series(axes[3, 7], time, trna_fraction_produced, 'RNA mass fraction produced',
				timeline, filtered, downsample=downsample)
			plt.tight_layout()
			return axes

		# Plot all time series data
		subplots(set(), downsample=1)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		# Downsample for less data and better illustrator load
		axes = subplots(set())
		exportFigure(plt, plotOutDir, f'{plotOutFileName}_downsampled', metadata)

		# Trim axes from all data for easier comparison across runs
		trim_axes(axes)
		exportFigure(plt, plotOutDir, f'{plotOutFileName}_trimmed', metadata)

		# Set axes limits for easier comparison across runs and filter time
		# points without all cells for smoother traces
		axes = subplots(filtered)
		trim_axes(axes)
		exportFigure(plt, plotOutDir, f'{plotOutFileName}_filtered', metadata)

		# Plot histograms of data
		_, axes = plt.subplots(4, 6, figsize=(20, 15))
		self.plot_hist(axes[0, 0], growth_rate, 0, 2, 'Growth rate\n(1/hr)')
		self.plot_hist(axes[1, 0], rna_growth, 0, 2, 'RNA growth rate\n(1/hr)')
		self.plot_hist(axes[2, 0], protein_growth, 0, 2, 'Protein growth rate\n(1/hr)')
		self.plot_hist(axes[3, 0], small_mol_growth, 0, 2, 'Small mol growth rate\n(1/hr)')
		self.plot_hist(axes[0, 1], rnap_elong_rate, 0, 100, 'RNAP elongation rate\n(nt/s)', sf=1)
		self.plot_hist(axes[1, 1], rnap_fraction_active, 0, 1, 'RNAP active fraction')
		self.plot_hist(axes[2, 1], ribosome_elong_rate, 0, 25, 'Ribosome elongation rate\n(AA/s)', sf=1)
		self.plot_hist(axes[3, 1], ribosome_fraction_active, 0, 1, 'Ribosome active fraction')
		self.plot_hist(axes[0, 2], fraction_charged[:, 0], 0, 1, f'Fraction charged\n{aa_ids[0][:-3]} tRNA')
		self.plot_hist(axes[1, 2], fraction_charged[:, 10], 0, 1, f'Fraction charged\n{aa_ids[10][:-3]} tRNA')
		self.plot_hist(axes[2, 2], aa_conc[:, 0], 0, 10, f'{aa_ids[0][:-3]} concentration\n(mM)')
		self.plot_hist(axes[3, 2], aa_conc[:, 10], 0, 1, f'{aa_ids[10][:-3]} concentration\n(mM)')
		self.plot_hist(axes[0, 3], ppgpp_conc, 0, 300, 'ppGpp concentration\n(uM)', n_bins=100, sf=1)
		self.plot_hist(axes[1, 3], fraction_charged, 0, 1, 'Fraction charged')
		self.plot_hist(axes[2, 3], aa_conc, 0, 100, 'Amino acid concentrations\n(mM)')  # TODO: handle log scale and variable ranges
		self.plot_hist(axes[3, 3], rna_fraction_prob, 0, 1, 'RNA fraction\nsynthesis probability')
		self.plot_hist(axes[0, 4], rp_ratio, 0, 1, 'RNA/protein\nmass fraction')
		self.plot_hist(axes[1, 4], rna_fraction, 0, 0.15, 'RNA mass fraction')
		self.plot_hist(axes[2, 4], protein_fraction, 0, 0.3, 'Protein mass fraction')
		self.plot_hist(axes[0, 5], rnap_conc, 0, 10, 'RNAP conc\n(uM)')
		self.plot_hist(axes[1, 5], rnap_output, 0, 0.2, 'RNAP output\n(mM NTPs/s)')
		self.plot_hist(axes[2, 5], ribosome_conc, 0, 40, 'Ribosome conc\n(uM)')
		self.plot_hist(axes[3, 5], ribosome_output, 0, 1, 'Ribosome output\n(mM AA/s)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, f'{plotOutFileName}_hist', metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
