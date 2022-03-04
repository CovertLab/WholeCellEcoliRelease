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
from models.ecoli.analysis.single.ribosome_limitation import calculate_ribosome_excesses
from models.ecoli.sim.variants import ppgpp_conc, ppgpp_limitations
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns, read_stacked_bulk_molecules
from wholecell.io.tablereader import TableReader
from wholecell.analysis.plotting_tools import COLORS_COLORBLIND as COLORS


MEAN = 'mean'
STD = 'std'


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def plot_data(self, axes, x, y, yerr, xlabel, ylabel, condition_labels, conditions,
			factors, condition_base_factor):
		raw_ax, norm_ax = axes
		for condition, color in zip(np.unique(conditions), COLORS):
			mask = conditions == condition
			raw_ax.errorbar(x[mask], y[mask], yerr=yerr[mask], fmt='o', color=color,
				label=condition_labels[condition])

			if condition_base_factor is not None and condition_base_factor.get(condition) in factors[mask]:
				ref_idx = factors[mask] == condition_base_factor[condition]
				ref_val = y[mask][ref_idx]
				norm_ax.errorbar(x[mask], y[mask] / ref_val, yerr=yerr[mask] / ref_val,
					fmt='o', color=color, label=condition_labels[condition])
				norm_ax.axhline(1, linestyle='--', color='k', linewidth=0.5)

		raw_ax.set_ylabel(ylabel, fontsize=8)
		norm_ax.set_ylabel(f'Normalized\n{ylabel}', fontsize=8)
		for ax in axes:
			ax.set_xlabel(xlabel, fontsize=8)
			ax.tick_params(labelsize=6)
			self.remove_border(ax)

	def plot_overlays(self, data, keys, x, twinx=False, normalize=4):
		_, (raw_aw, norm_ax) = plt.subplots(1, 2, figsize=(8, 2.5))
		create_axis = twinx
		for key, color in zip(keys, COLORS):
			y = data[key][MEAN]
			yerr = data[key][STD]
			raw_aw.errorbar(x, y, yerr=yerr, fmt='o', color=color, label=key)

			y_norm = y[normalize] if len(y) > normalize else 1
			norm_ax.errorbar(x, y / y_norm, yerr=yerr / y_norm, fmt='o', color=color, label=key)

			if twinx:
				raw_aw.set_ylabel(key)

			if create_axis:
				raw_aw = plt.twinx(raw_aw)
				create_axis = False

		norm_ax.set_ylabel('Normalized values')
		norm_ax.legend(fontsize=6)
		self.remove_border(raw_aw)
		self.remove_border(norm_ax)
		plt.tight_layout()

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

		variants = self.ap.get_variants()
		n_variants = len(variants)

		conditions = np.zeros(n_variants, int)
		factors = np.zeros(n_variants)
		data = {}
		if metadata.get('variant') == 'ppgpp_limitations':
			split_index = ppgpp_limitations.plot_split
		else:
			split_index = ppgpp_conc.split_index

		for i, variant in enumerate(variants):
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
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
			def add_data(key, time_series):
				entry = data.get(key, dict(mean=np.zeros(n_variants), std=np.zeros(n_variants)))
				entry[MEAN][i] = time_series.mean()
				entry[STD][i] = time_series.std()
				data[key] = entry

			add_data('ppgpp', ppgpp)
			add_data('growth_rate', growth_rate)
			add_data('rna_to_protein', rna_to_protein)
			add_data('elong_rate', elong_rate)
			add_data('ribosome_output', ribosome_output)
			add_data('aa_output', aa_output)
			add_data('ribosome_capacity', ribosome_capacity)
			add_data('aa_capacity', aa_capacity)
			add_data('fraction_charged', fraction_charged)
			add_data('ribosome_saturation', ribosome_saturation)
			add_data('aa_saturation', aa_saturation)
			add_data('aa_conc', aa_conc)
			add_data('excess_rna', excess[:, 0])
			add_data('excess_protein', excess[:, 1])
			add_data('rrna_synth_fraction', synth_fractions[:, 0])
			add_data('rprotein_synth_fraction', synth_fractions[:, 1])
			add_data('enzyme_synth_fraction', synth_fractions[:, 2])
			add_data('rprotein_protein_fraction', protein_fractions[:, 0])
			add_data('enzyme_protein_fraction', protein_fractions[:, 1])

		# Compile for plots
		if metadata.get('variant') == 'ppgpp_limitations':
			x = factors
			xlabel = 'Adjustment factor'
			condition_labels = {c: c for c in conditions}
			condition_base_factor = None
		else:
			x = data['ppgpp'][MEAN]
			xlabel = 'ppGpp conc (uM)'
			condition_labels = {i: c for i, c in enumerate(sim_data.ordered_conditions)}
			condition_base_factor = dict(zip(ppgpp_conc.CONDITIONS, ppgpp_conc.BASE_FACTOR))
		labels = {
			'growth_rate': 'Growth rate (1/hr)',
			'elong_rate': 'Elongation rate (AA/s)',
			'rna_to_protein': 'RNA/protein',
			'ribosome_output': 'Ribosome output (mM/s)',
			'aa_output': 'AA output (mM/s)',
			'ribosome_capacity': 'Ribosome capacity (mM/s)',
			'aa_capacity': 'AA capacity (mM/s)',
			'fraction_charged': 'Fraction charged',
			'ribosome_saturation': 'Ribosome saturation',
			'aa_saturation': 'AA synthesis average saturation',
			'aa_conc': 'AA conc (mM)',
			'excess_rna': 'rRNA excess',
			'excess_protein': 'rProtein excess',
			'rrna_synth_fraction': 'rRNA synth fraction\n(per rRNA, rProtein, enzymes mass)',
			'rprotein_synth_fraction': 'rProtein synth fraction\n(per rRNA, rProtein, enzymes mass)',
			'enzyme_synth_fraction': 'Enzyme synth fraction\n(per rRNA, rProtein, enzymes mass)',
			'rprotein_protein_fraction': 'rProtein fraction\n(per protein mass)',
			'enzyme_protein_fraction': 'Enzyme fraction\n(per protein mass)',
			}

		# Create plots
		n_subplots = len(labels)
		_, axes = plt.subplots(n_subplots, 2, figsize=(8, 2.5 * n_subplots))

		## Bar plots of cell properties
		for i, (key, ylabel) in enumerate(labels.items()):
			self.plot_data(axes[i, :], x, data[key][MEAN], data[key][STD],
				xlabel, ylabel, condition_labels, conditions, factors, condition_base_factor)

		## Formating for plots
		axes[0, 0].legend(fontsize=6)
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')

		# Plots specific to paper
		## Output
		keys = ['ribosome_output', 'aa_output']
		self.plot_overlays(data, keys, x)
		exportFigure(plt, plotOutDir, plotOutFileName + '_output', metadata)
		plt.close('all')

		## Capacity
		keys = ['ribosome_capacity', 'aa_capacity']
		self.plot_overlays(data, keys, x, twinx=True)
		exportFigure(plt, plotOutDir, plotOutFileName + '_capacity', metadata)
		plt.close('all')

		## Excess
		keys = ['excess_rna', 'aa_conc']
		self.plot_overlays(data, keys, x, twinx=True)
		exportFigure(plt, plotOutDir, plotOutFileName + '_excess', metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
