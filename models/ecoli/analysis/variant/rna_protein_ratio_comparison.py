"""
Plots the relationship between growth rate and RNA/protein ratios to compare
with experimental data in Scott et al., Science (2010).

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/29/19
"""

from __future__ import absolute_import, print_function, division

from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath
from wholecell.utils.sparkline import whitePadSparklineAxis

# First generation (counting from zero) from which to gather ratio values.
# If fewer generations were run, this script quits early without plotting
# anything.
FIRST_GENERATION = 2

FIGSIZE = (4, 4)

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

		# Get errorbar plot
		plt.figure(figsize=FIGSIZE)

		plt.style.use('seaborn-deep')
		color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
		marker_styles = ['o', '^', 'x']
		labels = ['basal', 'anaerobic', '+AA']

		ax = plt.subplot2grid((1, 1), (0, 0))

		for i in range(3):
			ax.errorbar(
				all_growth_rates[i].mean(),
				all_rna_to_protein_ratios[i].mean(),
				yerr=all_rna_to_protein_ratios[i].std(), color=color_cycle[0],
				mec=color_cycle[0], marker=marker_styles[i], markersize=8,
				mfc='white', linewidth=1, capsize=2, label=labels[i])

		# Add linear plot proposed in Scott et al. (2010)
		x_linear = np.linspace(0.05, 1.95, 100)
		y_linear = x_linear/4.5 + 0.087
		ax.plot(x_linear, y_linear, linewidth=2, color=color_cycle[2])

		ax.set_xlim([0, 2])
		ax.set_ylim([0, 0.7])
		ax.get_yaxis().get_major_formatter().set_useOffset(False)
		ax.get_xaxis().get_major_formatter().set_useOffset(False)

		whitePadSparklineAxis(ax)

		ax.tick_params(which='both', bottom=True, left=True,
			top=False, right=False, labelbottom=True, labelleft=True)

		ax.set_xlabel("Growth rate $\lambda$ (hour$^{-1}$)")
		ax.set_ylabel("RNA/protein mass ratio")
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		# Get clean version of errorbar plot
		ax.set_xlabel("")
		ax.set_ylabel("")
		ax.set_yticklabels([])
		ax.set_xticklabels([])
		exportFigure(plt, plotOutDir, plotOutFileName + "_clean", metadata)

		plt.close("all")

		# Get scatter version of plot
		plt.figure(figsize=FIGSIZE)
		ax = plt.subplot2grid((1, 1), (0, 0))

		options = {
			"edgecolors": color_cycle[0], "alpha": 0.25, "s": 20
			}

		ax.scatter(
			all_growth_rates[0], all_rna_to_protein_ratios[0],
			facecolors="none", marker="o", label=labels[0], **options)
		ax.scatter(
			all_growth_rates[1], all_rna_to_protein_ratios[1],
			facecolors="none", marker="^", label=labels[1], **options)
		ax.scatter(
			all_growth_rates[2], all_rna_to_protein_ratios[2],
			marker="x", label=labels[2], **options)

		x_linear = np.linspace(0.05, 2.45, 100)
		y_linear = x_linear/4.5 + 0.087
		ax.plot(x_linear, y_linear, linewidth=2, color=color_cycle[2])

		ax.set_xlim([0, 2.5])
		ax.set_ylim([0, 0.8])
		ax.get_yaxis().get_major_formatter().set_useOffset(False)
		ax.get_xaxis().get_major_formatter().set_useOffset(False)

		whitePadSparklineAxis(ax)

		ax.tick_params(which='both', bottom=True, left=True,
			top=False, right=False, labelbottom=True, labelleft=True)

		ax.set_xlabel("Growth rate $\lambda$ (hour$^{-1}$)")
		ax.set_ylabel("RNA/protein mass ratio")
		exportFigure(plt, plotOutDir, plotOutFileName + "_scatter", metadata)

if __name__ == "__main__":
	Plot().cli()
