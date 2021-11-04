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
	def plot_hist(self, ax, data, min_val, max_val, label, n_bins=100, sf=2):
		mean = data.mean()

		# Plot data
		patch = ax.hist(data, bins=n_bins, range=(min_val, max_val), alpha=0.7, histtype='step')[-1][0]
		color = patch.get_edgecolor()
		ax.axvline(mean, color=color, linestyle='--', linewidth=0.5)

		# Add mean +/- std text
		ax.text(mean, ax.get_ylim()[1], f'{mean:.{sf}f} +/- {data.std():.{sf+1}f}', fontsize=6, color=color)

		# Format axes
		self.remove_border(ax)
		ax.tick_params(labelsize=6)
		ax.set_xlabel(label, fontsize=8)
		ax.set_ylabel('Count', fontsize=8)

	def plot_sim(self, axes, input_dir):
		cell_paths = AnalysisPaths(input_dir, variant_plot=True).get_cells()
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
		synth_prob_per_cistron = read_stacked_columns(cell_paths, 'RnaSynthProb', 'rna_synth_prob_per_cistron',
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
		rrna_fraction_prob = synth_prob_per_cistron[:, cistron_data['is_rRNA']].sum(1)
		trna_fraction_prob = synth_prob_per_cistron[:, cistron_data['is_tRNA']].sum(1)
		mrna_fraction_prob = synth_prob_per_cistron[:, cistron_data['is_mRNA']].sum(1)

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
		self.plot_hist(axes[2, 2], aa_conc[:, 0], 0, 5, f'{aa_ids[0][:-3]} concentration\n(mM)')
		self.plot_hist(axes[3, 2], aa_conc[:, 10], 0, 1, f'{aa_ids[10][:-3]} concentration\n(mM)')
		self.plot_hist(axes[0, 3], ppgpp_conc, 0, 300, 'ppGpp concentration\n(uM)', sf=1)
		self.plot_hist(axes[1, 3], rrna_fraction_prob, 0, 1, 'rRNA fraction\nsynthesis probability')
		self.plot_hist(axes[2, 3], trna_fraction_prob, 0, 1, 'tRNA fraction\nsynthesis probability')
		self.plot_hist(axes[3, 3], mrna_fraction_prob, 0, 1, 'mRNA fraction\nsynthesis probability')

	def do_plot(self, inputDir1, plotOutDir, plotOutFileName, inputDir2, unused, metadata):
		_, axes = plt.subplots(4, 4, figsize=(15, 15))

		self.plot_sim(axes, inputDir1)
		self.plot_sim(axes, inputDir2)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
