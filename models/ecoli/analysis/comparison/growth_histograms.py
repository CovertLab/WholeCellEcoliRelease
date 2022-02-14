"""
Compare histograms for growth related properties.

TODO:
- Shares a lot of code with growth_time_series cohort plot (data load) - turn into a function?
"""

import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import comparisonAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	def plot_hist(self, ax, data, min_val, max_val, xlabel, label, n_bins=100, sf=2):
		mean = data.mean()

		# Plot data
		patch = ax.hist(data, bins=n_bins, range=(min_val, max_val), alpha=0.7, histtype='step',
			label=f'{label}: {mean:.{sf}f} +/- {data.std():.{sf+1}f}')[-1][0]
		color = patch.get_edgecolor()
		ax.axvline(mean, color=color, linestyle='--', linewidth=0.5)

		# Format axes
		self.remove_border(ax)
		ax.tick_params(labelsize=6)
		ax.set_xlabel(xlabel, fontsize=6)
		ax.set_ylabel('Count', fontsize=6)

	def plot_sim(self, axes, input_dir, sim_label):
		sim_data = self.read_sim_data_file(input_dir)

		transcription = sim_data.process.transcription
		ppgpp_id = sim_data.molecule_ids.ppGpp
		rnap_id = sim_data.molecule_ids.full_RNAP
		aa_ids = sim_data.molecule_groups.amino_acids
		ribosome_subunit_ids = [sim_data.molecule_ids.s30_full_complex, sim_data.molecule_ids.s50_full_complex]
		uncharged_trna_names = transcription.rna_data['id'][transcription.rna_data['is_tRNA']]
		charged_trna_names = transcription.charged_trna_names
		aa_from_trna = transcription.aa_from_trna.T
		cistron_data = transcription.cistron_data
		rna_fractions = ['is_rRNA', 'is_tRNA', 'is_mRNA']
		convert_to_fraction = lambda x: np.vstack([
			x[:, cistron_data[fraction]].sum(1)
			for fraction in rna_fractions
			]).T

		ap = AnalysisPaths(input_dir, variant_plot=True)
		for variant in ap.get_variants():
			cell_paths = ap.get_cells(variant=[variant])

			# Load attributes
			unique_molecule_reader = TableReader(os.path.join(cell_paths[0], 'simOut', 'UniqueMoleculeCounts'))
			unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
			rnap_idx = unique_molecule_ids.index('active_RNAP')
			ribosome_idx = unique_molecule_ids.index('active_ribosome')

			# Load data
			growth_function = lambda x: np.diff(x, axis=0) / x[:-1]
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
			rrna_fraction_prob, trna_fraction_prob, mrna_fraction_prob = read_stacked_columns(
				cell_paths, 'RnaSynthProb', 'rna_synth_prob_per_cistron',
				remove_first=True, ignore_exception=True, fun=convert_to_fraction).T
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

			label = f'Sim {sim_label}, var {variant}'
			self.plot_hist(axes[0, 0], growth_rate, 0, 2, 'Growth rate\n(1/hr)', label)
			self.plot_hist(axes[1, 0], rna_growth, 0, 2, 'RNA growth rate\n(1/hr)', label)
			self.plot_hist(axes[2, 0], protein_growth, 0, 2, 'Protein growth rate\n(1/hr)', label)
			self.plot_hist(axes[3, 0], small_mol_growth, 0, 2, 'Small mol growth rate\n(1/hr)', label)
			self.plot_hist(axes[0, 1], rnap_elong_rate, 0, 100, 'RNAP elongation rate\n(nt/s)', label, sf=1)
			self.plot_hist(axes[1, 1], rnap_fraction_active, 0, 1, 'RNAP active fraction', label)
			self.plot_hist(axes[2, 1], ribosome_elong_rate, 0, 25, 'Ribosome elongation rate\n(AA/s)', label, sf=1)
			self.plot_hist(axes[3, 1], ribosome_fraction_active, 0, 1, 'Ribosome active fraction', label)
			self.plot_hist(axes[0, 2], fraction_charged[:, 0], 0, 1, f'Fraction charged\n{aa_ids[0][:-3]} tRNA', label)
			self.plot_hist(axes[1, 2], fraction_charged[:, 10], 0, 1, f'Fraction charged\n{aa_ids[10][:-3]} tRNA', label)
			self.plot_hist(axes[2, 2], aa_conc[:, 0], 0, 10, f'{aa_ids[0][:-3]} concentration\n(mM)', label)
			self.plot_hist(axes[3, 2], aa_conc[:, 10], 0, 1, f'{aa_ids[10][:-3]} concentration\n(mM)', label)
			self.plot_hist(axes[0, 3], ppgpp_conc, 0, 300, 'ppGpp concentration\n(uM)', label, sf=1)
			self.plot_hist(axes[1, 3], rrna_fraction_prob, 0, 1, 'rRNA fraction\nsynthesis probability', label)
			self.plot_hist(axes[2, 3], trna_fraction_prob, 0, 1, 'tRNA fraction\nsynthesis probability', label)
			self.plot_hist(axes[3, 3], mrna_fraction_prob, 0, 1, 'mRNA fraction\nsynthesis probability', label)

	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		_, axes = plt.subplots(4, 4, figsize=(15, 15))

		self.plot_sim(axes, reference_sim_dir, 1)
		self.plot_sim(axes, input_sim_dir, 2)

		for ax in axes.flatten():
			ax.legend(fontsize=6)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
