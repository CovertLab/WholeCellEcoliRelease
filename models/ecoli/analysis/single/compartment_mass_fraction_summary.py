from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader

COLORS_256 = [ # From colorbrewer2.org, qualitative 8-class set 1
	[166,206,227],
	[31,120,180],
	[178,223,138],
	[255,222,0],
	[251,154,153],
	[227,26,28],
	[253,191,111],
	[255,127,0],
	[202,178,214],
	[106,61,154],
	[51,160,44]
	]

COLORS = [
	[colorValue/255. for colorValue in color]
	for color in COLORS_256
	]


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(
			self, simOutDir, plotOutDir, plotOutFileName,
			simDataFile, validationDataFile, metadata):
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		main_reader = TableReader(os.path.join(simOutDir, "Main"))

		cell = mass.readColumn("cellMass")

		projection = mass.readColumn("projection_mass")
		cytosol = mass.readColumn("cytosol_mass")
		extracellular = mass.readColumn("extracellular_mass")
		membrane = mass.readColumn("membrane_mass")
		outer_membrane = mass.readColumn("outer_membrane_mass")
		periplasm = mass.readColumn("periplasm_mass")
		pilus = mass.readColumn("pilus_mass")
		inner_membrane = mass.readColumn("inner_membrane_mass")
		flagellum = mass.readColumn("flagellum")

		initialTime = main_reader.readAttribute("initialTime")
		t = (main_reader.readColumn("time") - initialTime) / 60.

		masses = np.vstack([
			projection,
			cytosol,
			extracellular,
			membrane,
			outer_membrane,
			periplasm,
			pilus,
			inner_membrane,
			flagellum,
			]).T
		fractions = (masses / cell[:, None]).mean(axis=0)

		mass_labels = ["Projection", "Cytosol", "Extracellular", "Membrane",
			"Outer Membrane", "Periplasm", "Pilus", "Inner Membrane",
			"Flagellum"]
		legend = [
			'{} ({:.3e})'.format(label, fraction)
			for label, fraction in zip(mass_labels, fractions)
			] + ['Total cell mass']

		plt.figure(figsize = (8.5, 11))
		plt.gca().set_prop_cycle('color', COLORS)

		plt.plot(t, masses / masses[0, :], linewidth=2)
		plt.plot(t, cell / cell[0], color='k', linestyle=':')

		plt.title("Cell mass by compartment (average fraction of total cell mass in parentheses)")
		plt.xlabel("Time (min)")
		plt.ylabel("Mass (normalized by t = 0 min)")
		plt.legend(legend, loc="best")

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
