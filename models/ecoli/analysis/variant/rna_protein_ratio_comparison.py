"""
Plots the relationship between growth rate and RNA/protein ratios to compare
with experimental data in Scott et al., Science (2010).

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/29/19
"""

from __future__ import absolute_import, print_function, division

import cPickle
from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath, units

# First generation (counting from zero) from which to gather ratio values.
# If fewer generations were run, this script quits early without plotting
# anything.
FIRST_GENERATION = 2

FIGSIZE = (7.5, 7.5)

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata["variant"] != "condition":
			print('This analysis only runs for the "condition" variant.')
			return

		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		n_gens = ap.n_generation
		variants = ap.get_variants()

		if n_gens - 1 < FIRST_GENERATION:
			print('Not enough generations to plot.')
			return

		all_growth_rates = []
		all_rna_to_protein_ratios = []

		for variant in variants:
			doubling_times = np.zeros(0)
			variant_rna_to_protein_ratios = np.zeros(0)
			
			all_cells = ap.get_cells(
				variant=[variant],
				generation=range(FIRST_GENERATION, n_gens))

			if len(all_cells) == 0:
				continue

			for simDir in all_cells:
				try:
					simOutDir = os.path.join(simDir, "simOut")
					mass = TableReader(os.path.join(simOutDir, "Mass"))
					rna_mass = mass.readColumn("rnaMass")
					protein_mass = mass.readColumn("proteinMass")
					
					time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

					doubling_times = np.hstack(
						(doubling_times, (time[-1] - time[0])/3600.)
						)
					
					variant_rna_to_protein_ratios = np.hstack(
						(variant_rna_to_protein_ratios, rna_mass.mean()/protein_mass.mean())
						)
				except:
					continue

			variant_growth_rates = np.log(2)/doubling_times

			all_growth_rates.append(variant_growth_rates)
			all_rna_to_protein_ratios.append(variant_rna_to_protein_ratios)

		plt.figure(figsize=FIGSIZE)

		plt.style.use('seaborn-deep')
		color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

		for i in range(3):
			plt.errorbar(
				all_growth_rates[i].mean(),
				all_rna_to_protein_ratios[i].mean(),
				yerr=all_rna_to_protein_ratios[i].std(),
				color=color_cycle[0], marker='o', markersize=5, linewidth=1,
				capsize=2)

		# Add linear plot proposed in Scott et al. (2010)
		x_linear = np.linspace(0, 3, 100)
		y_linear = 0.23*x_linear + 0.09
		plt.plot(x_linear, y_linear, linewidth=2, color=color_cycle[2])

		plt.xlim([0, 3])
		plt.ylim([0, 1.6])
		plt.xlabel("Growth rate $\lambda$ (hour$^{-1}$)")
		plt.ylabel("RNA/protein mass ratio")
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
