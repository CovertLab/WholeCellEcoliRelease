"""
Compares tRNA charging kinetic parameter optimization solutions with
experimental measurements.
"""

import pickle
import os

from six.moves import cPickle

from matplotlib import pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants, units

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


Jakubowski_compared_measurements = {
	'GLNS-MONOMER[c]': 0.12,
	'ARGS-MONOMER[c]': 3.3,
	'GLURS-MONOMER[c]': 4,
	'ILES-MONOMER[c]': 0.3,
	'LEUS-MONOMER[c]': 1,
	'LYSS-CPLX[c]': 2.5,
	'PHES-CPLX[c]': 1.7,
	'THRS-CPLX[c]': 0.2,
	'VALS-MONOMER[c]': 6,
	}

color_measured = '#bdbdbd'
color_optimized = '#6b6ecf'
color_est_lower_limit = '#636363'
print_stats = True

class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = pickle.load(f)
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validation_data_file, 'rb') as f:
			validation_data = pickle.load(f)

		amino_acids = sim_data.molecule_groups.amino_acids
		i = amino_acids.index('L-SELENOCYSTEINE[c]')
		_ = amino_acids.pop(i)

		synthetases = [sim_data.relation.amino_acid_to_synthetase[amino_acid]
			for amino_acid in amino_acids]

		measurements = validation_data.trna_synthetase_kinetics
		kinetics = sim_data.relation.trna_charging_kinetics
		jakubowski = validation_data.jakubowski1984

		width = 0.6

		# Rank by most to least agreeing
		highest_report_by_others = []
		optimized = []
		for amino_acid in amino_acids:
			synthetase = sim_data.relation.amino_acid_to_synthetase[amino_acid]

			# Optimized by this study
			optimized_k_cat = sim_data.relation.synthetase_to_k_cat[synthetase]\
				.asNumber(1/units.s)
			optimized.append(optimized_k_cat)

			# Measurements
			others = []
			for substrate in measurements[synthetase].keys():
				others += [
					k.asNumber(1/units.s)
					for k in measurements[synthetase][substrate]['k_cat']]
			if synthetase in Jakubowski_compared_measurements:
				others.append(Jakubowski_compared_measurements[synthetase])

			# Estimates by Jakubowski & Goldman
			y = jakubowski[amino_acid]['Activity of synthetase in vivo']
			if isinstance(y, str):
				pass
			else:
				y = y.asNumber(1/units.s)
				others.append(y)

			# Highest report by others
			highest_report_by_others.append(max(others))

		# Rank synthetases by their agreement with others
		fold_change = np.array(optimized) / np.array(highest_report_by_others)
		rank = np.argsort(fold_change)
		synthetases = (np.array(synthetases)[rank]).tolist()
		boundary = np.where(fold_change[rank] > 10)[0][0] - 0.5

		fold_change = fold_change[rank]
		print('Fold change (optimized / highest other estimate or measurement)')
		i = synthetases.index('HISS-CPLX[c]')
		print(f'\tHisRS:\t{fold_change[i]:.2f}')
		i = synthetases.index('ARGS-MONOMER[c]')
		print(f'\tArgRS:\t{fold_change[i]:.2f}')

		# Curated in vitro measurements
		patches = []
		for x, synthetase in enumerate(synthetases):

			measured_k_cats = []
			for substrate in measurements[synthetase].keys():
				measured_k_cats += [k.asNumber(1/units.s)
					for k in measurements[synthetase][substrate]['k_cat']]

			if synthetase in Jakubowski_compared_measurements:
				measured_k_cats.append(Jakubowski_compared_measurements[synthetase])

			y = min(measured_k_cats)
			height = max(measured_k_cats) - y

			# To make single measurement viewable, give small height
			if height == 0:
				height = y / 30

			patches.append(
				Rectangle(
					(x - width / 2, y),
					width,
					height,
					)
				)
		curation_collection = PatchCollection(
			patches,
			facecolor=color_measured,
			edgecolor='none',
			alpha=0.8,
			)

		# Jakubowski and Goldman in vivo estimates
		x_values = []
		y_values = []
		for amino_acid in amino_acids:
			synthetase = sim_data.relation.amino_acid_to_synthetase[amino_acid]
			x = synthetases.index(synthetase)
			y = jakubowski[amino_acid]['Activity of synthetase in vivo']
			if isinstance(y, str):
				continue
			else:
				y = y.asNumber(1/units.s)
			x_values.append(x)
			y_values.append(y)

		# WCM Optimization
		patches = []
		for x, synthetase in enumerate(synthetases):

			k_cats = []
			for key in kinetics.keys():
				k_cats.append(
					kinetics[key]['synthetase_to_k_cat'][synthetase]
					.asNumber(1/units.s))

			y = min(k_cats)
			height = max(k_cats) - y
			patches.append(
				Rectangle(
					(x - width / 2, y),
					width,
					height,
					)
				)

		optimization_collection = PatchCollection(
			patches,
			facecolor=color_optimized,
			edgecolor='none',
			alpha=0.4,
			)

		x_values_wcm = []
		y_values_wcm = []
		for x, synthetase in enumerate(synthetases):
			x_values_wcm.append(x)
			y_values_wcm.append(
				sim_data.relation.synthetase_to_k_cat[synthetase]
				.asNumber(1/units.s))

		# Print statistics
		if print_stats:

			jakubowski_estimates = []
			optimized_estimates = []
			for amino_acid in amino_acids:
				synthetase = sim_data.relation.amino_acid_to_synthetase[amino_acid]

				# Estimates by Jakubowski & Goldman
				x = jakubowski[amino_acid]['Activity of synthetase in vivo']
				if isinstance(x, str):
					continue
				else:
					x = x.asNumber(1/units.s)
					jakubowski_estimates.append(x)

				# Optimized by this study
				optimized_estimates.append(sim_data.relation.synthetase_to_k_cat[synthetase].asNumber(1/units.s))

			print(f'Jakubowski estimates median:\t{np.median(jakubowski_estimates):.2f}')
			print(f'Optimized estimates median:\t{np.median(optimized_estimates):.2f}')

			highest_measurements = []
			optimized_estimates = []
			for amino_acid in amino_acids:
				synthetase = sim_data.relation.amino_acid_to_synthetase[amino_acid]

				# Measurements
				k_cats = []
				for substrate in measurements[synthetase].keys():
					k_cats += [k.asNumber(1/units.s)
						for k in measurements[synthetase][substrate]['k_cat']]
				if synthetase in Jakubowski_compared_measurements:
					k_cats.append(Jakubowski_compared_measurements[synthetase])
				highest_measurements.append(max(k_cats))

				# Optimized by this study
				optimized_estimates.append(sim_data.relation.synthetase_to_k_cat[synthetase].asNumber(1/units.s))

			fold_change = np.divide(optimized_estimates, highest_measurements)
			print(f'Average (median) fold change between optimized and highest measured: {np.median(fold_change):.2f}')


		# Plot
		fontsize = 9
		fig, ax = plt.subplots(1, 1, figsize=(6.5, 2.5))
		
		ax.add_collection(curation_collection)
		ax.add_collection(optimization_collection)
		ax.plot([], [], color='k', linewidth=0.5, label='in vitro measurement (refered by Jakubowski and Goldman, 1984)')
		ax.scatter(x_values, y_values, color=color_est_lower_limit, marker='^', label='in vivo estimate of lower limit activity (Jakubowski and Goldman, 1984)', s=6)
		ax.scatter(x_values_wcm, y_values_wcm, color=color_optimized, marker='o', label='Value used in this study', s=6)

		# Plot measurements Jakubowski and Goldman compared to
		for synthetase in Jakubowski_compared_measurements.keys():
			x = synthetases.index(synthetase)
			y = Jakubowski_compared_measurements[synthetase]
			ax.plot([x - width / 2, x + width / 2], [y, y], color='k', linewidth=0.5)

		# # Plot boundary (1 order of magnitude agreement)
		# ax.axvline(boundary, color='k', linewidth=0.5, linestyle='--')

		ax.set_xlabel('tRNA Synthetase', fontsize=fontsize)
		ax.set_ylabel('k_cat (1/s)', fontsize=fontsize)
		ax.set_yscale('log')
		ax.set_xticks(np.arange(len(synthetases)))
		ax.set_xticklabels(
			[synthetases_to_abbreviation[synthetase] for synthetase in synthetases],
			ha='center',
			va='top',
			rotation=90,
			fontsize=7,
			)
		ax.tick_params(axis='both', which='major', labelsize=7)

		# # Legend
		# handles, labels = ax.get_legend_handles_labels()
		# ax.legend(handles=[
		# 	mpatches.Patch(color=color_measured, alpha=0.4, label='Range of in vitro measurements'),
		# 	mpatches.Patch(color=color_optimized, alpha=0.4, label='Range of in vivo estimates'),
		# 	] + handles,
		# 	bbox_to_anchor=(0.5, 1),
		# 	loc='lower center',
		# 	ncol=3,
		# 	fontsize=5,
		# 	)

		ax.set_xlim([-1, 20])
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.margins(0.05)

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
