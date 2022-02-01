from __future__ import absolute_import, division, print_function

from matplotlib import pyplot as plt
import numpy as np
import os

from models.ecoli.analysis import cohortAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader, TableReaderError
from six.moves import range


class Plot(cohortAnalysisPlot.CohortAnalysisPlot):
	def do_plot(self, variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells in each seed

		max_cells_in_gen = 0
		for genIdx in range(self.ap.n_generation):
			n_cells = len(self.ap.get_cells(generation = [genIdx]))

			if n_cells > max_cells_in_gen:
				max_cells_in_gen = n_cells

		# noinspection PyTypeChecker
		fig, axesList = plt.subplots(self.ap.n_generation,
			sharey=True, sharex=True, figsize=(6, 3*self.ap.n_generation))

		all_final_masses = []

		for genIdx in range(self.ap.n_generation):
			gen_cells = self.ap.get_cells(generation = [genIdx])
			gen_final_masses = []

			for simDir in gen_cells:
				simOutDir = os.path.join(simDir, "simOut")
				mass_path = os.path.join(simOutDir, "Mass")

				try:
					# Listeners used
					mass_reader = TableReader(mass_path)

					# Load data
					cellMass = mass_reader.readColumn("cellMass")

					# Get final mass
					gen_final_masses.append(cellMass[-1] / 1000.)

				except (TableReaderError, EnvironmentError) as e:
					# Skip sims that were not able to complete division
					print("Couldn't read the Table {}; maybe the cell didn't finish division; skipping this sim: {!r}"
						.format(mass_path, e))
					continue

			all_final_masses.append(np.array(gen_final_masses))


		# Plot histograms of final masses
		if self.ap.n_generation == 1:
			axesList = [axesList]

		for idx, axes in enumerate(axesList):
			gen_final_masses = all_final_masses[idx]

			if max_cells_in_gen > 1:
				axes.hist(
					gen_final_masses,
					int(np.ceil(np.sqrt(gen_final_masses.size)))
					)
			else:
				axes.plot(gen_final_masses, 1, 'x')
				axes.set_ylim([0, 2])

			axes.axvline(gen_final_masses.mean(), color='k', linestyle='dashed', linewidth=2)
			axes.text(
				gen_final_masses.mean(), 1,
				"Mean = %.2f, Std = %.2f, n = %d" %
				(gen_final_masses.mean(), gen_final_masses.std(), len(gen_final_masses), ))
			axes.set_ylabel("Generation %d" % (idx, ))

		axesList[-1].set_xlabel("Final mass (pg)")

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
