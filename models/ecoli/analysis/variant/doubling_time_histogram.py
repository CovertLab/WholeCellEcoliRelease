import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.analysis.plotting_tools import COLORS_COLORBLIND as COLORS


FONT_SIZE=9
MAX_CELL_LENGTH = 180  # filter sims that reach the max time of 180 min


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def hist(self, ax, data, xlabel, bin_width=1., xlim=None):
		for variant, variant_data in data.items():
			color = COLORS[variant % len(COLORS)]
			bins = max(1, int(np.ceil((variant_data.max() - variant_data.min()) / bin_width)))
			mean = variant_data.mean()
			std = variant_data.std()
			ax.hist(variant_data, bins, color=color, alpha=0.5,
				label=f'Var {variant}: {mean:.1f} +/- {std:.2f}')
			ax.axvline(mean, color=color, linestyle='--', linewidth=1)

		if xlim:
			ax.set_xlim(xlim)
		self.remove_border(ax)
		ax.set_xlabel(xlabel, fontsize=FONT_SIZE)
		ax.tick_params(labelsize=FONT_SIZE)
		ax.legend()

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		doubling_times = {}
		growth_rates = {}

		def downsample(x):
			"""Average every n_downsample points to one value to smooth and downsample"""
			n_downsample = 100
			if (extra_points := x.shape[0] % n_downsample) != 0:
				x = x[:-extra_points]
			return x.reshape(-1, n_downsample).mean(1).reshape(-1, 1)

		for variant in self.ap.get_variants():
			all_cells = self.ap.get_cells(variant=[variant], only_successful=True)
			if len(all_cells) == 0:
				continue

			dt = read_stacked_columns(all_cells, 'Main', 'time',
				fun=lambda x: (x[-1] - x[0]) / 60.).squeeze()
			doubling_times[variant] = dt[dt < MAX_CELL_LENGTH]
			growth_rates[variant] = read_stacked_columns(all_cells[dt < MAX_CELL_LENGTH], 'Mass', 'instantaneous_growth_rate',
				remove_first=True, fun=downsample).squeeze() * 3600.

		_, axes = plt.subplots(2, 1, figsize=(10, 10))

		self.hist(axes[0], doubling_times, 'Doubling Time (min)')
		self.hist(axes[1], growth_rates, 'Growth rates (1/hr)', bin_width=0.05)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		axes[0].set_xlim([15, 90])
		axes[1].set_xlim([0, 2.5])
		exportFigure(plt, plotOutDir, plotOutFileName + '_trimmed', metadata)


if __name__ == "__main__":
	Plot().cli()
