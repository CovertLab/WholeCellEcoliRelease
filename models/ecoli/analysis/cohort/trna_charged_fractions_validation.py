"""
Plots trna charged fractions.
"""

import pickle
import os

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.plotting_tools import trna_to_label, amino_acid_to_label

probe_to_trnas = {
	'ArgUCU': ['argU-tRNA[c]'],
	'ArgACG': ['argV-tRNA[c]', 'argY-tRNA[c]', 'argQ-tRNA[c]', 'argZ-tRNA[c]'],
	'ArgCCG': ['argX-tRNA[c]'],
	'ArgCCU': ['argW-tRNA[c]'],

	'HisGUG': ['hisR-tRNA[c]'],
	'LysUUU': ['RNA0-303[c]', 'RNA0-301[c]', 'lysV-tRNA[c]', 'lysW-tRNA[c]', 'RNA0-302[c]', 'lysT-tRNA[c]'],
	'AspGUC': ['aspV-tRNA[c]', 'aspU-tRNA[c]', 'aspT-tRNA[c]'],
	'GluUUC': ['gltU-tRNA[c]', 'gltV-tRNA[c]', 'gltT-tRNA[c]', 'gltW-tRNA[c]'],
	'AsnGUU': ['asnU-tRNA[c]', 'asnV-tRNA[c]', 'RNA0-304[c]', 'asnT-tRNA[c]'],
	'CysGCA': ['cysT-tRNA[c]'],
	'GlnU/CUG': ['glnX-tRNA[c]', 'glnW-tRNA[c]', 'glnV-tRNA[c]', 'glnU-tRNA[c]'],

	'SerCGA': ['serU-tRNA[c]'],
	'SerGCU': ['serV-tRNA[c]'],
	'SerGGA': ['serW-tRNA[c]', 'serX-tRNA[c]'],
	'SerUGA': ['serT-tRNA[c]'],
	
	'ThrGGU': ['thrT-tRNA[c]', 'thrV-tRNA[c]'],
	'ThrCGU': ['thrW-tRNA[c]'],
	'ThrUGU': ['thrU-tRNA[c]'],

	'AlaG/UGC': ['alaX-tRNA[c]', 'alaV-tRNA[c]', 'alaW-tRNA[c]', 'alaT-tRNA[c]', 'alaU-tRNA[c]'],
	
	'GlyGCC': ['glyY-tRNA[c]', 'glyX-tRNA[c]', 'glyV-tRNA[c]', 'glyW-tRNA[c]'],
	'GlyUCC': ['glyT-tRNA[c]'],
	'GlyCCC': ['glyU-tRNA[c]'],
	
	'IleGAU': ['ileT-tRNA[c]', 'ileU-tRNA[c]', 'ileV-tRNA[c]'],
	'IleCAU': ['RNA0-305[c]', 'ileX-tRNA[c]'],
	
	'LeuCAA': ['leuX-tRNA[c]'],
	'LeuCAG': ['leuQ-tRNA[c]', 'leuV-tRNA[c]', 'leuP-tRNA[c]', 'leuT-tRNA[c]'],
	'LeuUAA': ['leuZ-tRNA[c]'],
	'LeuUAG': ['leuW-tRNA[c]'],

	'MetCAUf': ['RNA0-306[c]', 'metY-tRNA[c]', 'metZ-tRNA[c]', 'metW-tRNA[c]'],
	'MetCAUm': ['metT-tRNA[c]', 'metU-tRNA[c]'],
	
	'PheGAA': ['pheU-tRNA[c]', 'pheV-tRNA[c]'],
	
	'ProUGG': ['proM-tRNA[c]'],
	'ProCGG': ['proK-tRNA[c]'],
	'ProGGG': ['proL-tRNA[c]'],
	
	'TrpCCA': ['trpT-tRNA[c]'],
	'TyrGUA': ['tyrT-tRNA[c]', 'tyrV-tRNA[c]', 'tyrU-tRNA[c]'],
	
	'ValGAC': ['valW-tRNA[c]', 'valV-tRNA[c]'],
	'ValUAC': ['RNA0-300[c]', 'valY-tRNA[c]', 'valU-tRNA[c]', 'valX-tRNA[c]', 'valT-tRNA[c]']}

amino_acid_abv_to_id = {
	'Ala': 'L-ALPHA-ALANINE[c]',
	'Arg': 'ARG[c]',
	'Asn': 'ASN[c]',
	'Asp': 'L-ASPARTATE[c]',
	'Cys': 'CYS[c]',
	'Glu': 'GLT[c]', # Glutamic acid: 
	'Gln': 'GLN[c]', # Glutamine: 
	'Gly': 'GLY[c]',
	'His': 'HIS[c]',
	'Ile': 'ILE[c]',
	'Leu': 'LEU[c]',
	'Lys': 'LYS[c]',
	'Met': 'MET[c]',
	'Phe': 'PHE[c]',
	'Pro': 'PRO[c]',
	'Ser': 'SER[c]',
	'Thr': 'THR[c]',
	'Trp': 'TRP[c]',
	'Tyr': 'TYR[c]',
	'Val': 'VAL[c]',
}

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells()

		# Get molecules
		rna_data = sim_data.process.transcription.rna_data
		raw_trnas = rna_data['id'][rna_data['is_tRNA']].tolist()
		charged_trnas = sim_data.process.transcription.charged_trna_names

		# Remove selC-tRNA
		selC = raw_trnas.index('selC-tRNA[c]')
		_ = raw_trnas.pop(selC)
		_ = charged_trnas.pop(selC)

		# Order alphabetically for plotting
		trnas = [trna_to_label(trna) for trna in raw_trnas]
		indexes = np.argsort(trnas)
		trnas = np.array(trnas)[indexes]
		raw_trnas = np.array(raw_trnas)[indexes]
		charged_trnas = np.array(charged_trnas)[indexes]

		# Get number of molecules
		(n_charged_trnas, n_free_trnas) = read_stacked_bulk_molecules(
			cell_paths, (charged_trnas, raw_trnas), remove_first=True)

		# Calculate charged fraction for each tRNA
		f_charged = n_charged_trnas / (n_free_trnas + n_charged_trnas) * 100

		# Calculate charged fraction for each amino acid family
		f_charged_families = []
		for amino_acid in sim_data.molecule_groups.amino_acids:
			if amino_acid in ['CYS[c]', 'HIS[c]', 'TRP[c]', 'L-SELENOCYSTEINE[c]']:
				continue
			_trna_ids = sim_data.relation.amino_acid_to_trnas[amino_acid]
			_trna_indexes = [np.where(raw_trnas == trna)[0][0] for trna in _trna_ids]
			_trna_labels = [trna_to_label(trna) for trna in _trna_ids]

			_n_free_trnas_sum = n_free_trnas[:, _trna_indexes].sum(axis=1)
			_n_charged_trnas_sum = n_charged_trnas[:, _trna_indexes].sum(axis=1)
			_f_charged = (_n_charged_trnas_sum
				/ (_n_free_trnas_sum + _n_charged_trnas_sum)
				* 100)
			f_charged_families.append(_f_charged)

		################################################################
		# Plot
		fig, axes = plt.subplots(2, 1, figsize=(11, 4.5),
			gridspec_kw={'height_ratios': [9, 1]})

		# Add spaces between amino acid famililes
		xticklabels = []
		previous = 'ala'
		for trna in trnas:
			amino_acid = trna[:3]

			if amino_acid != previous:
				if previous not in ['cys', 'his', 'trp']:
					xticklabels.append('All')
				xticklabels.append('')
				previous = amino_acid

			xticklabels.append(trna)
		xticklabels.append('All')
		xticklabels.append('')
		positions = np.where(np.array(xticklabels) != '')[0]
		positions_trnas = np.where(np.logical_and(
			np.array(xticklabels) != 'All',
			np.array(xticklabels) != '',
			))[0]

		# Formatting
		fontsize = 7
		widths = 0.8
		alpha = 0.7

		# Plot simulated fraction charged
		def format_violins(violins, color):
			for pc in violins['bodies']:
				pc.set_facecolor(color)
				pc.set_linewidths(0.5)
				pc.set_alpha(alpha)

			for key in ['cbars', 'cmaxes', 'cmins']:
				violins[key].set_color(color)
				violins[key].set_linewidth(0.5)
			return

		ax = axes[0]
		color = '#6b6ecf'
		violins = ax.violinplot(
			f_charged, positions=positions_trnas, widths=widths)
		format_violins(violins, color)

		# Plot simulated fraction charged for each amino acid family
		positions_total = np.where(np.array(xticklabels) == 'All')[0]
		violins = ax.violinplot(
			f_charged_families, positions=positions_total, widths=widths)
		format_violins(violins, color)

		handle_simulated = mpatches.Patch(
			facecolor=color, linewidth=0, alpha=alpha,
			label='Simulated')

		################################################################
		# Plot measurements, Sorensen 2001
		sorensen_measurements = [
			# tRNA2ARG 74 +/- 1.6% charged (argV)
			['argV', 74, 1.6],
			# tRNA1THR 71.5 +/- 0.36% charged (thrV)
			['thrV', 71.5, 0.36],
			# tRNA1LEU 77% charged (leuP, leuQ, leuT, leuV)
			['leuP', 77, 0],
			['leuQ', 77, 0],
			['leuT', 77, 0],
			['leuV', 77, 0],
			# tRNAHIS 50% or 70 - 80% charged (hisR)
			['hisR', 50, 0],
			['hisR', 75, 5],
			]

		color = '#77cecc' # blue
		for measurement in sorensen_measurements:
			trna, y, yerr = measurement
			x = xticklabels.index(trna)
			ax.errorbar(x, y, yerr=yerr,
				marker='s', mfc=color, mec=color, ms=3,
				ecolor=color, elinewidth=0.5, mew=0.5, capsize=1.2)

		handle_sorensen = Line2D([0], [0],
			marker='s', mfc=color, mec=color, mew=0.5, ms=3, linewidth=0,
			label='Sorensen 2001')

		################################################################
		# Plot measurements, Kruger and Sorensen 1998
		kruger_measurements = [
			# tRNA2Glu 70.2 +/- 2.6% (gltT, gltU, gltV, gltW)
			['gltT', 70.2, 2.6],
			['gltU', 70.2, 2.6],
			['gltV', 70.2, 2.6],
			['gltW', 70.2, 2.6],
			# tRNALys 75.8 +/- 2.3% (lysT, lysV, lysW, lysY, lysZ, lysQ)
			['lysT', 75.8, 2.3],
			['lysV', 75.8, 2.3],
			['lysW', 75.8, 2.3],
			['lysY', 75.8, 2.3],
			['lysZ', 75.8, 2.3],
			['lysQ', 75.8, 2.3],
			# tRNA1Gln 69.3 +/- 1.8% (glnU, glnW)
			['glnU', 69.3, 1.8],
			['glnW', 69.3, 1.8],
			]

		color = '#b5cf6b' # green
		for measurement in kruger_measurements:
			trna, y, yerr = measurement
			x = xticklabels.index(trna)
			ax.errorbar(x, y, yerr=yerr,
				marker='s', mfc=color, mec=color, ms=3,
				ecolor=color, elinewidth=0.5, mew=0.5, capsize=1.2)

		handle_kruger = Line2D([0], [0],
			marker='s', mfc=color, mec=color, mew=0.5, ms=3, linewidth=0,
			label='Kruger et al. 1998')

		################################################################
		# Plot measurements, Avcilar-Kucukgoze et al. 2016
		avcilar_kucukgoze_measurements = [
			['ArgUCU', 1.311, 0.320],
			['ArgACG', 0.409, 0.264],
			['ArgCCG', 1.180, 0.180],
			['ArgCCU', 0.508, 0.273],
			['HisGUG', 0.362, 0.026],
			['LysUUU', 0.591, 0.008],
			['AspGUC', 0.756, 0.031],
			['GluUUC', 1.285, 0.518],
			['AsnGUU', 0.438, 0.001],
			['CysGCA', 0.347, 0.115],
			['GlnU/CUG', 1.097, 0.087],
			['SerCGA', 0.347, 0.077],
			['SerGCU', 0.541, 0.242],
			['SerGGA', 0.523, 0.171],
			['SerUGA', 0.401, 0.123],
			['ThrGGU', 0.476, 0.181],
			['ThrCGU', 0.339, 0.109],
			['ThrUGU', 0.489, 0.040],
			['AlaG/UGC', 0.545, 0.511],
			['GlyGCC', 0.595, 0.050],
			['GlyUCC', 0.533, 0.148],
			['GlyCCC', 0.686, 0.101],
			['IleGAU', 0.223, 0.086],
			['IleCAU', 0.561, 0.097],
			['LeuCAA', 0.454, 0.134],
			['LeuCAG', 0.880, 0.026],
			['LeuUAA', 0.463, 0.049],
			['LeuUAG', 0.503, 0.084],
			['MetCAUf', 0.545, 0.150],
			['MetCAUm', 0.879, 0.213],
			['PheGAA', 0.321, 0.210],
			['ProUGG', 0.753, 0.123],
			['ProCGG', 0.888, 0.064],
			['ProGGG', 0.721, 0.193],
			['TrpCCA', 0.778, 0.018],
			['TyrGUA', 0.599, 0.034],
			['ValGAC', 0.846, 0.219],
			['ValUAC', 0.795, 0.150],
			]

		color = '#e7ba52' # yellow
		for measurement in avcilar_kucukgoze_measurements:
			probe, y, yerr = measurement

			_raw_trnas = probe_to_trnas[probe]
			_trnas = [trna_to_label(trna) for trna in _raw_trnas]
			amino_acid = amino_acid_abv_to_id[probe[:3]]
			amino_acid_all_trnas = (
				sim_data.relation.amino_acid_to_trnas[amino_acid])
			k_M_Ts = [sim_data.relation.trna_to_K_T[trna]
				.asNumber(units.umol/units.L) for trna in _raw_trnas]

			# Measurements of single tRNAs
			if len(_trnas) == 1:
				x = xticklabels.index(_trnas[0])

			# Measurements of all tRNAs in the same amino acid family
			elif len(_trnas) == len(amino_acid_all_trnas):
				_trna_indexes = [xticklabels.index(trna) for trna in _trnas]
				x = max(_trna_indexes) + 1

			# Measurements of tRNAs with the same K_M_T (and thus the
			# same f charged in our model)
			elif len(set(k_M_Ts)) == 1:
				x = [xticklabels.index(trna) for trna in _trnas]
				y = np.ones_like(x) * y

			else:
				print(measurement, _trnas)
				continue

			ax.errorbar(x, y * 100, yerr=yerr * 100,
				marker='s', mfc=color, mec=color, ms=3, ls='none',
				ecolor=color, elinewidth=0.5, mew=0.5, capsize=1.2)

		handle_avcilar_kucukgoze = Line2D([0], [0],
			marker='s', mfc=color, mec=color, mew=0.5, ms=3, linewidth=0,
			label='Avcilar-Kucukgoze et al. 2016')		

		################################################################
		# Plot measurements, Yegian et al. 1966
		yegian_measurements = [
			['L-ALPHA-ALANINE[c]', 83],
			['L-ASPARTATE[c]', 39],
			['GLT[c]', 90],
			['GLY[c]', 75],
			['ILE[c]', 44],
			['LEU[c]', 97],
			['MET[c]', 102],
			['PHE[c]', 82],
			['PRO[c]', 73],
			['SER[c]', 33],
			['THR[c]', 80],
			['TYR[c]', 56],
			['VAL[c]', 97],
			['ARG[c]', 87],
			['HIS[c]', 61],
			['LYS[c]', 91],
			]

		color = '#d6616b' # red
		for measurement in yegian_measurements:
			amino_acid, y = measurement
			_trnas = [trna_to_label(trna) for trna
				in sim_data.relation.amino_acid_to_trnas[amino_acid]]
			_trna_indexes = [xticklabels.index(trna) for trna in _trnas]

			if amino_acid in ['CYS[c]', 'HIS[c]', 'TRP[c]']:
				x = _trna_indexes[0]
			else:
				x = max(_trna_indexes) + 1
			ax.errorbar(x, y,
				marker='s', mfc=color, mec=color, ms=3,
				ecolor=color, elinewidth=0.5, mew=0.5, capsize=1.2)

		handle_yegian = Line2D([0], [0],
			marker='s', mfc=color, mec=color, mew=0.5, ms=3, linewidth=0,
			label='Yegian et al. 1966')

		################################################################
		# Formatting
		ax.axhline(100, color='k', linewidth=0.5, linestyle='--')
		ax.legend(handles=[
			handle_simulated,
			handle_sorensen,
			handle_kruger,
			handle_avcilar_kucukgoze,
			handle_yegian,
			],
			loc='upper right', ncol=2, fontsize=fontsize)

		ax.set_xticks(positions)
		ax.set_xticklabels(
			[label for label in xticklabels if label != ''], rotation=90)
		ax.set_ylabel('Aminoacylated Percent', fontsize=fontsize)
		ax.tick_params(axis='both', which='major', labelsize=fontsize)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# Unify x axes
		for ax in axes:
			ax.set_xlim([-1, max(positions) + 1])

		################################################################
		# Delineate amino acid families
		ax = axes[1]
		for spine in ax.spines.values():
			spine.set_visible(False)
		ax.set_xticks([])
		ax.set_yticks([])

		space = 0.4
		for amino_acid in sim_data.molecule_groups.amino_acids:
			if amino_acid == 'L-SELENOCYSTEINE[c]':
				continue

			_trnas = [trna_to_label(trna)
				for trna in sim_data.relation.amino_acid_to_trnas[amino_acid]]
			x = [xticklabels.index(trna) for trna in _trnas]
			x_min = min(x) - space
			x_max = max(x) + space

			if amino_acid not in ['CYS[c]', 'HIS[c]', 'TRP[c]']:
				x_max += 1
			ax.plot([x_min, x_max], [1, 1], color='k', linewidth=1)

			if len(_trnas) <= 3 or amino_acid == 'GLN[c]':
				text = _trnas[0][0].upper() + _trnas[0][1:3]
			else:
				text = amino_acid_to_label[amino_acid]
			ax.text(np.mean(x), 0.95, text, fontsize=fontsize, ha='center',
				va='top')

		plt.tight_layout(h_pad=-0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
