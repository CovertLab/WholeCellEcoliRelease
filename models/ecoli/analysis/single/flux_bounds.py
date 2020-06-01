"""
Plot upper and lower flux targets

@author: Sophie Landon
@organization: University of Bristol
@date: Created 14/5/2020
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		mainListener = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = mainListener.readAttribute("initialTime")
		time = mainListener.readColumn("time") - initialTime

		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		targetFluxesUpper = enzymeKineticsReader.readColumn('targetFluxesUpper')
		targetFluxesLower = enzymeKineticsReader.readColumn('targetFluxesLower')
		actualFluxes = enzymeKineticsReader.readColumn('actualFluxes')
		reactionNames = [x[:18] for x in enzymeKineticsReader.readAttribute('constrainedReactions')]


		fig = plt.figure(figsize = (34, 34))
		ax_full = fig.add_subplot(111)
		ax_full.set_xlabel('Time', fontsize = 48)
		ax_full.set_ylabel('Flux', fontsize = 48)
		ax_full.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

		n_fluxes = actualFluxes.shape[1]
		n_rows = int(np.ceil(np.sqrt(n_fluxes)))
		n_cols = int(np.ceil(n_fluxes / n_rows))

		for idx in range(n_fluxes):

			ax = fig.add_subplot(n_rows, n_cols, idx + 1)
			plt.title(reactionNames[idx], fontsize = 7)

			plt.plot(time / 60., targetFluxesUpper[:, idx], linewidth = 2, label = 'Upper Target Flux')
			plt.plot(time / 60., targetFluxesLower[:, idx], linewidth = 2, label = 'Lower Target Flux')
			plt.plot(time / 60., actualFluxes[:, idx], linewidth = 2, label = 'Actual Flux')
			plt.tick_params(labelsize=0, top=False, bottom=False, right=False, left=False)

		handles, labels = ax.get_legend_handles_labels()
		ax_full.legend(handles, labels, bbox_to_anchor = (0.35, 1.15), fontsize = 48)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
