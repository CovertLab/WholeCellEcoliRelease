"""
Plots tRNA synthetase stability.
"""
from six.moves import cPickle

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.fireworks.firetasks.parca import ParcaTask
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


synthetases_to_abbreviation = {
	'ALAS-CPLX[c]': 'AlaRS',
	'ARGS-MONOMER[c]': 'ArgRS',
	'ASNS-CPLX[c]': 'AsnRS',
	'ASPS-CPLX[c]': 'AspRS',
	'CYSS-MONOMER[c]': 'CysRS',
	'GLURS-MONOMER[c]': 'GluRS',
	'GLNS-MONOMER[c]': 'GlnRS',
	'GLYS-CPLX[c]': 'GlyRS',
	'HISS-CPLX[c]': 'HisRS',
	'ILES-MONOMER[c]': 'IleRS',
	'LEUS-MONOMER[c]': 'LeuRS',
	'LYSS-CPLX[c]': 'LysRS',
	'METG-CPLX[c]': 'MetRS',
	'PHES-CPLX[c]': 'PheRS',
	'PROS-CPLX[c]': 'ProRS',
	'SERS-CPLX[c]': 'SerRS',
	'THRS-CPLX[c]': 'ThrRS',
	'TRPS-CPLX[c]': 'TrpRS',
	'TYRS-CPLX[c]': 'TyrRS',
	'VALS-MONOMER[c]': 'ValRS',
	}

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		# Molecules
		synthetases = []
		for amino_acid in sim_data.molecule_groups.amino_acids:
			if amino_acid == 'L-SELENOCYSTEINE[c]':
				continue
			synthetases.append(sim_data.relation.amino_acid_to_synthetase[amino_acid])

		# Simulated tRNA synthetase concentrations
		fig, axes = plt.subplots(2, 1, figsize=(6, 6.5))
		color = '#6b6ecf'
		color_accent = '#636363'
		color_ticks = '#969696'

		titles = [
			'M9 Minimal Media +0.4% Glucose (Aerobic)',
			'M9 Minimal Media +0.4% Glucose + Amino Acids (Aerobic)',
			]

		positions = np.arange(len(synthetases))
		widths = 0.8

		n_increments = 5
		sweep_levels = np.arange(2, n_increments + 1)

		for variant, ax, title in zip(variants, axes, titles):

			# Load modified variant sim_data
			with open(ap.get_variant_kb(variant), 'rb') as f:
				variant_sim_data = pickle.load(f)
			condition = variant_sim_data.condition

			# Get cells
			cell_paths = ap.get_cells(variant=[variant])

			# Load data
			(n_synthetases,) = read_stacked_bulk_molecules(
				cell_paths, (synthetases,))
			volume = (units.L
				* 1e-15
				* read_stacked_columns(cell_paths, 'Mass', 'cellVolume'))
			simulated = (1
				/ sim_data.constants.n_avogadro
				/ volume
				* n_synthetases
				).asNumber(units.umol/units.L)

			# Plot
			ax.set_title(title, fontsize=11)
			violins = ax.violinplot(
				simulated,
				positions=positions,
				widths=widths,
				)

			for pc in violins['bodies']:
				pc.set_facecolor(color)
				pc.set_edgecolor('none')
				pc.set_alpha(0.5)

			for key in ['cbars', 'cmaxes', 'cmins']:
				violins[key].set_color(color)
				violins[key].set_linewidth(0.5)

			# Indicate sweep levels
			minimum = simulated.min(axis=0)

			for i in range(len(synthetases)):
				x = positions[i]
				for sweep_level in sweep_levels:
					y = minimum[i] / n_increments * sweep_level
					ax.plot([x - widths/4, x + widths/4], [y, y], linewidth=0.5, color=color_ticks)
					ax.text(x + widths/3, y, sweep_level, fontsize=2, va='center', ha='left')
				ax.plot([x, x], [0, minimum[i]], linewidth=0.5, color=color_ticks)

		# Format and save plot
		axes[1].set_xlabel('tRNA Synthetases', fontsize=9)
		xticklabels = [synthetases_to_abbreviation[synthetase]
			for synthetase in synthetases]
		for ax in axes:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)

			ax.set_ylabel('Concentration (uM)', fontsize=9)
			ax.set_xticks(range(len(synthetases)))
			ax.set_xticklabels(xticklabels, rotation=90, fontsize=7)
			ax.tick_params(axis='both', which='major', labelsize=7)
			ax.axhline(0, color='k', linestyle='--', linewidth=0.5)

		plt.tight_layout(h_pad=2)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
