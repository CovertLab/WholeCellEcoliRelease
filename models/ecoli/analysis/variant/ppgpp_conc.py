"""
Compare cell properties at varying levels of ppGpp concentration.  Useful with
with the ppgpp_conc variant and in comparison with data presented in Zhu et al.
2019. https://academic.oup.com/nar/article/47/9/4684/5420536.
"""

import os
import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.analysis.single.ribosome_limitation import calculate_ribosome_excesses
from models.ecoli.sim.variants.ppgpp_conc import BASE_FACTOR, CONDITIONS, split_index
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, read_stacked_bulk_molecules
from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import COLORS_SMALL


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def plot_data(self, axes, ppgpp, y, yerr, ylabel, condition_labels, conditions, factors):
		condition_base_factor = dict(zip(CONDITIONS, BASE_FACTOR))

		raw_ax, norm_ax = axes
		for condition, color in zip(np.unique(conditions), COLORS_SMALL):
			mask = conditions == condition
			raw_ax.errorbar(ppgpp[mask], y[mask], yerr=yerr[mask], fmt='o', color=color,
				label=condition_labels[condition])

			if condition_base_factor[condition] in factors[mask]:
				ref_idx = factors[mask] == condition_base_factor[condition]
				ref_val = y[mask][ref_idx]
				norm_ax.errorbar(ppgpp[mask], y[mask] / ref_val, yerr=yerr[mask] / ref_val,
					fmt='o', color=color, label=condition_labels[condition])
				norm_ax.axhline(1, linestyle='--', color='k', linewidth=0.5)

		raw_ax.set_ylabel(ylabel, fontsize=8)
		norm_ax.set_ylabel(f'Normalized\n{ylabel}', fontsize=8)
		for ax in axes:
			ax.set_xlabel('ppGpp conc', fontsize=8)
			ax.tick_params(labelsize=6)
			self.remove_border(ax)

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		metabolism = sim_data.process.metabolism
		transcription = sim_data.process.transcription
		aa_enzyme_ids = metabolism.aa_enzymes
		get_enzymes = metabolism.get_pathway_enzyme_counts_per_aa
		fwd_kcats = metabolism.aa_kcats_fwd
		ribosome_subunit_ids = [sim_data.molecule_ids.s30_full_complex, sim_data.molecule_ids.s50_full_complex]
		aa_ids = sim_data.molecule_groups.amino_acids
		max_elong_rate = sim_data.process.translation.basal_elongation_rate
		uncharged_trna_names = transcription.rna_data['id'][transcription.rna_data['is_tRNA']]
		charged_trna_names = transcription.charged_trna_names
		aa_from_trna = transcription.aa_from_trna.T

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		conditions = np.zeros(n_variants, int)
		factors = np.zeros(n_variants)
		ppgpp_mean = np.zeros(n_variants)
		ppgpp_std = np.zeros(n_variants)
		growth_rate_mean = np.zeros(n_variants)
		growth_rate_std = np.zeros(n_variants)
		rna_to_protein_mean = np.zeros(n_variants)
		rna_to_protein_std = np.zeros(n_variants)
		elong_rate_mean = np.zeros(n_variants)
		elong_rate_std = np.zeros(n_variants)
		ribosome_output_mean = np.zeros(n_variants)
		ribosome_output_std = np.zeros(n_variants)
		aa_output_mean = np.zeros(n_variants)
		aa_output_std = np.zeros(n_variants)
		ribosome_capacity_mean = np.zeros(n_variants)
		ribosome_capacity_std = np.zeros(n_variants)
		aa_capacity_mean = np.zeros(n_variants)
		aa_capacity_std = np.zeros(n_variants)
		fraction_charged_mean = np.zeros(n_variants)
		fraction_charged_std = np.zeros(n_variants)
		ribosome_saturation_mean = np.zeros(n_variants)
		ribosome_saturation_std = np.zeros(n_variants)
		aa_saturation_mean = np.zeros(n_variants)
		aa_saturation_std = np.zeros(n_variants)
		aa_conc_mean = np.zeros(n_variants)
		aa_conc_std = np.zeros(n_variants)
		excess_rna_mean = np.zeros(n_variants)
		excess_rna_std = np.zeros(n_variants)
		excess_protein_mean = np.zeros(n_variants)
		excess_protein_std = np.zeros(n_variants)
		rrna_synth_fraction_mean = np.zeros(n_variants)
		rrna_synth_fraction_std = np.zeros(n_variants)
		rprotein_synth_fraction_mean = np.zeros(n_variants)
		rprotein_synth_fraction_std = np.zeros(n_variants)
		enzyme_synth_fraction_mean = np.zeros(n_variants)
		enzyme_synth_fraction_std = np.zeros(n_variants)
		rprotein_protein_fraction_mean = np.zeros(n_variants)
		rprotein_protein_fraction_std = np.zeros(n_variants)
		enzyme_protein_fraction_mean = np.zeros(n_variants)
		enzyme_protein_fraction_std = np.zeros(n_variants)
		for i, variant in enumerate(variants):
			all_cells = ap.get_cells(variant=[variant], only_successful=True)
			conditions[i], factors[i] = split_index(variant)

			unique_molecule_reader = TableReader(os.path.join(all_cells[0], 'simOut', 'UniqueMoleculeCounts'))
			unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
			ribosome_idx = unique_molecule_ids.index('active_ribosome')

			# Read data from listeners
			ppgpp = read_stacked_columns(all_cells, 'GrowthLimits', 'ppgpp_conc')
			elong_rate = read_stacked_columns(all_cells, 'RibosomeData', 'effectiveElongationRate')
			growth_rate = read_stacked_columns(all_cells, 'Mass', 'instantaneous_growth_rate', remove_first=True) * 3600
			rna_mass = read_stacked_columns(all_cells, 'Mass', 'rnaMass')
			protein_mass = read_stacked_columns(all_cells, 'Mass', 'proteinMass')
			aas_elongated = read_stacked_columns(all_cells, 'GrowthLimits', 'aasUsed', remove_first=True)
			aas_supplied = read_stacked_columns(all_cells, 'GrowthLimits', 'aa_supply', remove_first=True)
			aa_saturation = read_stacked_columns(all_cells, 'GrowthLimits', 'aa_supply_fraction_fwd', remove_first=True).mean(1)
			time_step = read_stacked_columns(all_cells, 'Main', 'timeStepSec', remove_first=True).squeeze()
			counts_to_molar = read_stacked_columns(all_cells, 'EnzymeKinetics', 'countsToMolar', remove_first=True).squeeze()
			unique_mol_counts = read_stacked_columns(all_cells, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts', remove_first=True)
			enzymes, ribosome_subunits, aas, uncharged_trna_counts, charged_trna_counts = read_stacked_bulk_molecules(
				all_cells, (aa_enzyme_ids, ribosome_subunit_ids, aa_ids, uncharged_trna_names, charged_trna_names), remove_first=True)
			excess, synth_fractions, protein_fractions = calculate_ribosome_excesses(sim_data, all_cells)

			rna_to_protein = (rna_mass / protein_mass)
			ribosome_output = counts_to_molar * aas_elongated.sum(1) / time_step
			aa_output = counts_to_molar * aas_supplied.sum(1) / time_step
			aa_conc = counts_to_molar * aas.sum(1)

			fwd_enzymes, rev_enzymes = get_enzymes(enzymes)
			aa_capacity = fwd_enzymes @ fwd_kcats * counts_to_molar
			active_ribosome_counts = unique_mol_counts[:, ribosome_idx]
			inactive_ribosome_counts = ribosome_subunits.min(1)
			ribosome_capacity = (active_ribosome_counts + inactive_ribosome_counts) * counts_to_molar * max_elong_rate
			ribosome_saturation = elong_rate / max_elong_rate

			charged_trna_counts = charged_trna_counts @ aa_from_trna
			uncharged_trna_counts = uncharged_trna_counts @ aa_from_trna
			fraction_charged = charged_trna_counts / (uncharged_trna_counts + charged_trna_counts)

			# Calculate mean and std for each value
			ppgpp_mean[i] = ppgpp.mean()
			ppgpp_std[i] = ppgpp.std()
			growth_rate_mean[i] = growth_rate.mean()
			growth_rate_std[i] = growth_rate.std()
			rna_to_protein_mean[i] = rna_to_protein.mean()
			rna_to_protein_std[i] = rna_to_protein.std()
			elong_rate_mean[i] = elong_rate.mean()
			elong_rate_std[i] = elong_rate.std()
			ribosome_output_mean[i] = ribosome_output.mean()
			ribosome_output_std[i] = ribosome_output.std()
			aa_output_mean[i] = aa_output.mean()
			aa_output_std[i] = aa_output.std()
			ribosome_capacity_mean[i] = ribosome_capacity.mean()
			ribosome_capacity_std[i] = ribosome_capacity.std()
			aa_capacity_mean[i] = aa_capacity.mean()
			aa_capacity_std[i] = aa_capacity.std()
			fraction_charged_mean[i] = fraction_charged.mean()
			fraction_charged_std[i] = fraction_charged.std()
			ribosome_saturation_mean[i] = ribosome_saturation.mean()
			ribosome_saturation_std[i] = ribosome_saturation.std()
			aa_saturation_mean[i] = aa_saturation.mean()
			aa_saturation_std[i] = aa_saturation.std()
			aa_conc_mean[i] = aa_conc.mean()
			aa_conc_std[i] = aa_conc.std()
			excess_rna_mean[i] = excess[:, 0].mean()
			excess_rna_std[i] = excess[:, 0].std()
			excess_protein_mean[i] = excess[:, 1].mean()
			excess_protein_std[i] = excess[:, 1].std()
			rrna_synth_fraction_mean[i] = synth_fractions[:, 0].mean()
			rrna_synth_fraction_std[i] = synth_fractions[:, 0].std()
			rprotein_synth_fraction_mean[i] = synth_fractions[:, 1].mean()
			rprotein_synth_fraction_std[i] = synth_fractions[:, 1].std()
			enzyme_synth_fraction_mean[i] = synth_fractions[:, 2].mean()
			enzyme_synth_fraction_std[i] = synth_fractions[:, 2].std()
			rprotein_protein_fraction_mean[i] = protein_fractions[:, 0].mean()
			rprotein_protein_fraction_std[i] = protein_fractions[:, 0].std()
			enzyme_protein_fraction_mean[i] = protein_fractions[:, 1].mean()
			enzyme_protein_fraction_std[i] = protein_fractions[:, 1].std()

		condition_labels = sim_data.ordered_conditions

		# Create plots
		_, axes = plt.subplots(18, 2, figsize=(8, 45))

		## Bar plots of cell properties
		self.plot_data(axes[0, :], ppgpp_mean, growth_rate_mean, growth_rate_std,
			'Growth rate (1/hr)', condition_labels, conditions, factors)
		self.plot_data(axes[1, :], ppgpp_mean, elong_rate_mean, elong_rate_std,
			'Elongation rate (AA/s)', condition_labels, conditions, factors)
		self.plot_data(axes[2, :], ppgpp_mean, rna_to_protein_mean, rna_to_protein_std,
			'RNA/protein', condition_labels, conditions, factors)
		self.plot_data(axes[3, :], ppgpp_mean, ribosome_output_mean, ribosome_output_std,
			'Ribosome output (uM/s)', condition_labels, conditions, factors)
		self.plot_data(axes[4, :], ppgpp_mean, aa_output_mean, aa_output_std,
			'AA output (uM/s)', condition_labels, conditions, factors)
		self.plot_data(axes[5, :], ppgpp_mean, ribosome_capacity_mean, ribosome_capacity_std,
			'Ribosome capacity (uM/s)', condition_labels, conditions, factors)
		self.plot_data(axes[6, :], ppgpp_mean, aa_capacity_mean, aa_capacity_std,
			'AA capacity (uM/s)', condition_labels, conditions, factors)
		self.plot_data(axes[7, :], ppgpp_mean, fraction_charged_mean, fraction_charged_std,
			'Fraction charged', condition_labels, conditions, factors)
		self.plot_data(axes[8, :], ppgpp_mean, ribosome_saturation_mean, ribosome_saturation_std,
			'Ribosome saturation', condition_labels, conditions, factors)
		self.plot_data(axes[9, :], ppgpp_mean, aa_saturation_mean, aa_saturation_std,
			'AA synthesis average saturation', condition_labels, conditions, factors)
		self.plot_data(axes[10, :], ppgpp_mean, aa_conc_mean, aa_conc_std,
			'AA conc (uM)', condition_labels, conditions, factors)
		self.plot_data(axes[11, :], ppgpp_mean, excess_rna_mean, excess_rna_std,
			'rRNA excess', condition_labels, conditions, factors)
		self.plot_data(axes[12, :], ppgpp_mean, excess_protein_mean, excess_protein_std,
			'rProtein excess', condition_labels, conditions, factors)
		self.plot_data(axes[13, :], ppgpp_mean, rrna_synth_fraction_mean, rrna_synth_fraction_std,
			'rRNA synth fraction\n(per rRNA, rProtein, enzymes mass)', condition_labels, conditions, factors)
		self.plot_data(axes[14, :], ppgpp_mean, rprotein_synth_fraction_mean, rprotein_synth_fraction_std,
			'rProtein synth fraction\n(per rRNA, rProtein, enzymes mass)', condition_labels, conditions, factors)
		self.plot_data(axes[15, :], ppgpp_mean, enzyme_synth_fraction_mean, enzyme_synth_fraction_std,
			'Enzyme synth fraction\n(per rRNA, rProtein, enzymes mass)', condition_labels, conditions, factors)
		self.plot_data(axes[16, :], ppgpp_mean, rprotein_synth_fraction_mean, rprotein_synth_fraction_std,
			'rProtein fraction\n(per protein mass)', condition_labels, conditions, factors)
		self.plot_data(axes[17, :], ppgpp_mean, enzyme_synth_fraction_mean, enzyme_synth_fraction_std,
			'Enzyme fraction\n(per protein mass)', condition_labels, conditions, factors)

		axes[0, 0].legend(fontsize=6)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
