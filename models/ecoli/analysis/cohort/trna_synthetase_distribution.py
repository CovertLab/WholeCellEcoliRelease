"""
Plots abundance and distribution of tRNA synthetases.
"""

import pickle
import os

from scipy.stats import gamma
from scipy.special import gamma as gamma_function
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np

from models.ecoli.analysis import cohortAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import (exportFigure,
	read_bulk_molecule_counts, read_stacked_bulk_molecules, read_stacked_columns)
from wholecell.io.tablereader import TableReader

# Taniguchi et al. 2010
measurements = {
	'ILES-MONOMER[c]': {'mean': 64.061, 'sd': 22.707, 'a': 8.1, 'b': 8.148,},
	'CYSS-MONOMER[c]': {'mean': 80, 'sd': 30.499, 'a': 7.13, 'b': 12.072,},
	'LEUS-MONOMER[c]': {'mean': 96.516, 'sd': 32.569, 'a': 8.995, 'b': 11.149,},
	'GLNS-MONOMER[c]': {'mean': 67.294, 'sd': 21.511, 'a': 9.965, 'b': 6.959,},
	'SERS-CPLX[c]': {'mean': 168.148, 'sd': 66.683, 'a': 6.56, 'b': 27.179,},
	'THRS-CPLX[c]': {'mean': 783.577, 'sd': 320.924, 'a': 5.92, 'b': 132.368,},
	'ASPS-CPLX[c]': {'mean': 182.294, 'sd': 58.538, 'a': 9.78, 'b': 18.8,},
	'METG-CPLX[c]': {'mean': 79.581, 'sd': 29.678, 'a': 7.22, 'b': 11.174,},
	'HISS-CPLX[c]': {'mean': 75.695, 'sd': 24.312, 'a': 9.995, 'b': 7.959,},
	'ALAS-CPLX[c]': {'mean': 86.147, 'sd': 32.126, 'a': 7.875, 'b': 12.091,},
	'LYSS-CPLX[c]': {'mean': 118.608, 'sd': 47.894, 'a': 6.165, 'b': 19.351,},
	'TRPS-CPLX[c]': {'mean': 125.146, 'sd': 39.192, 'a': 10.29, 'b': 12.293,},
	'VALS-MONOMER[c]': {'mean': 103.936, 'sd': 30.611, 'a': 11.965, 'b': 9.178,},
}

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

def get_fit(x, k, theta):
	numerator = np.power(x, k-1) * np.exp(-x / theta)
	denominator = gamma_function(k) * (theta**k)
	return numerator / denominator

class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		ap = AnalysisPaths(variantDir, cohort_plot=True)
		cell_paths = ap.get_cells()

		# Load data
		amino_acids = sim_data.molecule_groups.amino_acids
		_ = amino_acids.pop(amino_acids.index('L-SELENOCYSTEINE[c]'))
		synthetases = [sim_data.relation.amino_acid_to_synthetase[amino_acid]
			for amino_acid in amino_acids]
		(n_synthetases,) = read_stacked_bulk_molecules(cell_paths, (synthetases,))

		# Plot
		fig, axes = plt.subplots(5, 4, figsize=(8, 8))
		axes = axes.reshape(-1)

		color_bar = '#9c9ede'
		color_measurement = '#e7ba52'

		for ax, synthetase in zip(axes, synthetases):

			# Get data
			i = synthetases.index(synthetase)
			data = n_synthetases[:, i]

			# Histogram
			bars = ax.hist(data)
			ax.cla()

			bins = bars[1]
			bin_size = bins[1] - bins[0]
			weights = np.ones(n_synthetases.shape[0]) / n_synthetases.shape[0] / bin_size

			x_min = np.floor(n_synthetases[:, i].min() / bin_size) * bin_size
			x_max = np.ceil(n_synthetases[:, i].max() / bin_size) * bin_size
			bins = np.arange(x_min, x_max, bin_size)
			bars = ax.hist(data, bins=bins, weights=weights, color=color_bar, edgecolor='w')
			xlim = ax.get_xlim()

			shape, loc, scale = gamma.fit(data, floc=0)

			# Validation gamma distribution (scaled to WCM mean)
			if synthetase in measurements:

				# Align validation mean to model mean
				factor = data.mean() / measurements[synthetase]['mean']
				mean = measurements[synthetase]['mean'] * factor
				std = measurements[synthetase]['sd'] * factor

				# Compute k and theta
				theta = std**2 / mean
				k = mean / theta

				# Plot gamma distribution
				x = np.linspace(x_min, x_max, num=100)
				y = get_fit(x, k, theta)
				ax.plot(x, y, color=color_measurement)

			# Format
			ax.set_title(synthetases_to_abbreviation[synthetase], fontsize=9)
			ax.set_xlabel('Number of Molecules', fontsize=7)
			ax.set_ylabel('Probability Density', fontsize=7)

			ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
			ax.yaxis.offsetText.set_fontsize(7)
			ax.tick_params(axis='both', which='major', labelsize=7)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)

		simulated = mpatches.Patch(color=color_bar, label='Simulated')
		measured = mlines.Line2D([], [], color=color_measurement, label='Scaled Measurement')
		# plt.suptitle('Probability Density Functions')
		fig.legend(handles=[simulated, measured], loc='upper right', fontsize=7)

		plt.tight_layout(rect=(0, 0, 1, 0.95))
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
