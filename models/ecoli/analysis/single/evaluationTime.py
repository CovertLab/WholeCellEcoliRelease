"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/10/2014
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


def subplot(gs, x, y, title, labels, sort=False):
	ax = plt.subplot(gs)
	if sort:
		idx = np.argsort(y[-1, :])[::-1]
	else:
		idx = np.arange(y.shape[1])

	# Determine legends from total time
	total_time = (y / 60).sum(axis=0)
	legend_labels = np.array(['{} ({:.2f})'.format(name, t)
		for name, t in zip(labels, total_time)])

	# Plot
	ax.semilogy(x, 1000 * y[:, idx])

	# Formatting
	ax.grid(True, which='major')
	ax.set_xlabel('Simulation time (min)')
	ax.set_ylabel('Evaluation time (ms)')
	ax.set_title(title)
	ax.legend(legend_labels[idx], bbox_to_anchor=(1,1), prop={'size':6})


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		evaluationTime = TableReader(os.path.join(simOutDir, "EvaluationTime"))
		mainReader = TableReader(os.path.join(simOutDir, "Main"))

		stateNames = evaluationTime.readAttribute("stateNames")
		processNames = evaluationTime.readAttribute("processNames")

		update_queries = evaluationTime.readColumn("updateQueries_times")
		partition = evaluationTime.readColumn("partition_times")
		merge = evaluationTime.readColumn("merge_times")
		calculate_request = evaluationTime.readColumn("calculateRequest_times")
		evolve_state = evaluationTime.readColumn("evolveState_times")

		initialTime = mainReader.readAttribute("initialTime")
		time = (mainReader.readColumn("time") - initialTime) / 60  # min

		plt.figure(figsize=(10, 10))
		gs = GridSpec(3, 2)
		subplot(gs[0, 0], time, update_queries, 'State.updateQueries', stateNames)
		subplot(gs[1, 0], time, partition, 'State.partition', stateNames)
		subplot(gs[2, 0], time, merge, 'State.merge', stateNames)
		subplot(gs[0, 1], time, calculate_request, 'Process.calculateRequest', processNames, sort=True)
		subplot(gs[1, 1], time, evolve_state, 'Process.evolveState', processNames, sort=True)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
