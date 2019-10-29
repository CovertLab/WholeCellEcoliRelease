"""
Gets scatter plots of birth volumes vs added volumes in different nutrient
conditions, to compare with experimental data in Wallden et al.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/19/19
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
from wholecell.utils.sparkline import whitePadSparklineAxis

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

		initial_volumes = []
		added_volumes = []

		for variant in variants:
			with open(ap.get_variant_kb(variant), 'rb') as f:
				sim_data = cPickle.load(f)

			cell_density = sim_data.constants.cellDensity

			initial_masses = np.zeros(0)
			final_masses = np.zeros(0)

			all_cells = ap.get_cells(variant=[variant], generation=gens)

			if len(all_cells) == 0:
				continue

			for simDir in all_cells:
				try:
					simOutDir = os.path.join(simDir, "simOut")
					mass = TableReader(os.path.join(simOutDir, "Mass"))
					cellMass = mass.readColumn("cellMass")

					initial_masses = np.hstack((initial_masses, cellMass[0]))
					final_masses = np.hstack((final_masses, cellMass[-1]))
				except:
					continue

			added_masses = final_masses - initial_masses

			initial_volume = initial_masses/cell_density.asNumber(units.fg/units.um**3)
			added_volume = added_masses/cell_density.asNumber(units.fg/units.um**3)

			initial_volumes.append(initial_volume)
			added_volumes.append(added_volume)

		plt.style.use('seaborn-deep')
		color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

		plt.figure(figsize=(4, 4))
		ax = plt.subplot2grid((1, 1), (0, 0))

		options = {
			"edgecolors": color_cycle[0], "alpha": 0.2, "s": 50, "clip_on": False
			}
		labels = ["minimal", "anaerobic", "minimal + AA"]

		ax.scatter(initial_volumes[2], added_volumes[2],
			marker="x", label=labels[2], **options)
		ax.scatter(initial_volumes[0], added_volumes[0],
			facecolors="none", marker="o", label=labels[0], **options)
		ax.scatter(initial_volumes[1], added_volumes[1],
			facecolors="none", marker="^", label=labels[1], **options)

		ax.set_xlim([0, 4])
		ax.set_ylim([0, 4])
		ax.set_xlabel("Birth Volume ($\mu m^3$)")
		ax.set_ylabel("Added Volume ($\mu m^3$)")
		ax.legend()

		ax.get_yaxis().get_major_formatter().set_useOffset(False)
		ax.get_xaxis().get_major_formatter().set_useOffset(False)

		whitePadSparklineAxis(ax)

		ax.tick_params(which='both', bottom=True, left=True,
			top=False, right=False, labelbottom=True, labelleft=True)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		# Get clean version of plot
		ax.set_xlabel("")
		ax.set_ylabel("")
		ax.set_yticklabels([])
		ax.set_xticklabels([])
		exportFigure(plt, plotOutDir, plotOutFileName + "_clean", metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
