"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		evaluationTime = TableReader(os.path.join(simOutDir, "EvaluationTime"))

		stateNames = evaluationTime.readAttribute("stateNames")
		processNames = evaluationTime.readAttribute("processNames")

		updateQueries_times = evaluationTime.readColumn("updateQueries_times")
		partition_times = evaluationTime.readColumn("partition_times")
		merge_times = evaluationTime.readColumn("merge_times")
		calculateRequest_times = evaluationTime.readColumn("calculateRequest_times")
		evolveState_times = evaluationTime.readColumn("evolveState_times")

		evaluationTime.close()

		mainReader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = mainReader.readAttribute("initialTime")
		time = mainReader.readColumn("time") - initialTime
		mainReader.close()

		plt.figure(figsize = (8.5, 11))

		plt.subplot(3, 2, 1)

		plt.semilogy(time / 60, updateQueries_times * 1000)
		plt.grid(True, which = "major")
		plt.xlabel("Simulation time (min)")
		plt.ylabel("Evaluation time (ms)")
		plt.title("State.updateQueries")
		plt.legend(stateNames, loc="best", bbox_to_anchor=(1,1), prop={'size':6})

		plt.subplot(3, 2, 2)

		idx = np.argsort(calculateRequest_times[-1, :])[::-1]
		plt.semilogy(time / 60, calculateRequest_times[:, idx] * 1000)
		plt.grid(True, which = "major")
		plt.xlabel("Simulation time (min)")
		plt.ylabel("Evaluation time (ms)")
		plt.title("Process.calculateRequest")
		plt.legend(np.array(processNames)[idx], loc="best", bbox_to_anchor=(1,1), prop={'size':6})

		plt.subplot(3, 2, 3)

		plt.semilogy(time / 60, partition_times * 1000)
		plt.grid(True, which = "major")
		plt.xlabel("Simulation time (min)")
		plt.ylabel("Evaluation time (ms)")
		plt.title("State.partition")
		plt.legend(stateNames, loc="best", bbox_to_anchor=(1,1), prop={'size':6})

		plt.subplot(3, 2, 4)

		idx = np.argsort(evolveState_times[-1, :])[::-1]
		plt.semilogy(time / 60, evolveState_times[:, idx] * 1000)
		plt.grid(True, which = "major")
		plt.xlabel("Simulation time (min)")
		plt.ylabel("Evaluation time (ms)")
		plt.title("Process.evolveState")
		plt.legend(np.array(processNames)[idx], loc="best", bbox_to_anchor=(1,1), prop={'size':6})

		plt.subplot(3, 2, 5)

		plt.semilogy(time / 60, merge_times * 1000)
		plt.grid(True, which = "major")
		plt.xlabel("Simulation time (min)")
		plt.ylabel("Evaluation time (ms)")
		plt.title("State.merge")
		plt.legend(stateNames, loc="best", bbox_to_anchor=(1,1), prop={'size':6})

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
