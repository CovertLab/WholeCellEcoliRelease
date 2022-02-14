from __future__ import annotations

import os

from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.analysis import comparisonAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants

COLORS = [
	'tab:cyan', 'tab:blue',
	'tab:orange', 'tab:brown',
	'tab:pink', 'tab:purple',
	'tab:olive', 'tab:green',
	'tab:gray', 'tab:red']


class Plot(comparisonAnalysisPlot.ComparisonAnalysisPlot):
	mass_names = ['dryMass', 'proteinMass', 'rRnaMass', 'mRnaMass', 'dnaMass']
	mass_labels = [name.replace('Mass', ' (fg)') for name in mass_names]

	def do_plot(self, reference_sim_dir, plotOutDir, plotOutFileName, input_sim_dir, unused, metadata):
		# noinspection PyTypeChecker
		_, self.axs = plt.subplots(len(self.mass_names), sharex=True)

		self.max_time = 0

		self.plot_wcm(reference_sim_dir, 0)
		self.plot_wcm(input_sim_dir, 1)

		for idx, ax in enumerate(self.axs):
			ax.set_xlim(0, self.max_time)
			ax.set_yticks(list(ax.get_ylim()))
			ax.set_ylabel(self.mass_labels[idx])

		self.axs[0].set_title("Cell mass fraction comparisons")
		plt.legend(bbox_to_anchor=(.8, 5), loc=2, borderaxespad=0., prop={'size': 6})
		self.axs[-1].set_xlabel("Time (hr)")
		plt.subplots_adjust(hspace=0.2, wspace=0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

	def plot_wcm(self, inputDir, num):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		all_cells = ap.get_cells()
		variant_was_labeled = set()  # label each variant once

		for sim_path in all_cells:
			variant_id = ap.get_cell_variant(sim_path)
			color = COLORS[(variant_id * 2 + num) % len(COLORS)]
			label = '_'  # a no-op label

			if variant_id not in variant_was_labeled:
				variant_was_labeled.add(variant_id)
				variant_dir = os.path.join(sim_path, '..', '..', '..')
				with open(os.path.join(variant_dir, constants.METADATA_DIR,
									   'short_name')) as f:
					label = next(f).strip()

			simOutDir = os.path.join(sim_path, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
			self.max_time = max(self.max_time, time[-1] / 3600)
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			for mass_idx, mass_name in enumerate(self.mass_names):
				mass_data = mass.readColumn(mass_name)
				self.axs[mass_idx].plot(
					time / 3600,
					mass_data,
					linestyle='solid' if num == 0 else 'dashed',
					linewidth=3 if num == 0 else 1,
					color=color,
					label=label)


if __name__ == "__main__":
	Plot().cli()
