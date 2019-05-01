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

		plt.figure(figsize=(5, 5))
		plt.scatter(initial_volumes[0], added_volumes[0], s=3, label="minimal")
		plt.scatter(initial_volumes[1], added_volumes[1], s=3, label="anaerobic")
		plt.scatter(initial_volumes[2], added_volumes[2], s=3, label="+AA")
		plt.xlim([0, 4])
		plt.ylim([0, 4])
		plt.xlabel("Birth Volume ($\mu m^3$)")
		plt.ylabel("Added Volume ($\mu m^3$)")
		plt.legend()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
