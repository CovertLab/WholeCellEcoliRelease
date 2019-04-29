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


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata["variant"] != "condition":
			print("This plot only runs for the 'condition' variant.")
			return

		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		gens = [2, 3]

		growth_rates = []
		rna_to_protein_ratios = []

		for variant in variants:
			doubling_times = np.zeros(0)
			rna_to_protein_ratio = np.zeros(0)
			
			all_cells = ap.get_cells(variant=[variant], generation=gens)

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
					
					rna_to_protein_ratio = np.hstack(
						(rna_to_protein_ratio, rna_mass.mean()/protein_mass.mean())
						)
				except:
					continue

			growth_rate = np.log(2)/doubling_times

			growth_rates.append(growth_rate)
			rna_to_protein_ratios.append(rna_to_protein_ratio)

		plt.figure(figsize=(7.5, 7.5))
		plt.scatter(growth_rates[0], rna_to_protein_ratios[0], s=3, color='c', label="minimal")
		plt.scatter(growth_rates[1], rna_to_protein_ratios[1], s=3, color='m', label="anaerobic")
		plt.scatter(growth_rates[2], rna_to_protein_ratios[2], s=3, color='k', label="+AA")
		plt.xlim([0, 3])
		plt.ylim([0, 0.75])
		plt.xlabel("Growth rate $\lambda$ (hour$^{-1}$)")
		plt.ylabel("RNA/protein mass ratio")
		plt.legend()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
