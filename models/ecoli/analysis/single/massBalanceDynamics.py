from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

THRESHOLD = 1e-5


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		mass = TableReader(os.path.join(simOutDir, "Mass"))

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = (
			TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
			) / 60.

		relProcessMassDifferences = mass.readColumn("relProcessMassDifferences")

		relProcessMassDifferences[relProcessMassDifferences == 0] = np.nan

		processNames = mass.readAttribute("processNames")

		mass.close()

		plt.figure(figsize = (8.5, 11))

		n_processes = len(processNames)

		n_cols = int(np.sqrt(n_processes))
		n_rows = int(np.ceil(n_processes/n_cols))

		axis = [time.min(), time.max(), np.nanmin(np.abs(relProcessMassDifferences)), np.nanmax(np.abs(relProcessMassDifferences))]

		for i, processName in enumerate(processNames):
			plt.subplot(n_rows, n_cols, i+1)

			series = relProcessMassDifferences[:, i].copy()
			t = time.copy()

			t = t[series != 0]
			series = series[series != 0]

			if np.any(series > 0):
				colors = np.zeros((series.size, 3), np.float64)

				colors[np.abs(series) > THRESHOLD, :] = (1.0, 0.0, 0.0)

				# markers = ["^" if value > 0 else "v" for value in series]

				plt.scatter(
					t,
					np.abs(series),
					color = colors,
					marker = '.'
					# marker = markers
					)

			plt.hlines(THRESHOLD,  axis[0], axis[1], "r", "dashed")

			plt.title(processName)

			# plt.xlabel("time (min)")
			# plt.ylabel("mass diff ($|(m_f - m_i)/(m_i)|$)")

			plt.yscale("log")

			plt.axis(axis)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
