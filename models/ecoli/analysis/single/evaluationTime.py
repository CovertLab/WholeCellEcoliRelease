from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import zip


def subplot(gs, x, y, title, labels, sort=False):
	assert y.ndim == 2, 'y.ndim={}, title={}, labels={}'.format(y.ndim, title, labels)
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
	ax.set_title(title + ' (total {:.2f} mins)'.format(total_time.sum()))
	ax.legend(legend_labels[idx], bbox_to_anchor=(1,1), prop={'size':6},
		loc='upper left')


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		evaluationTime = TableReader(os.path.join(simOutDir, "EvaluationTime"))
		mainReader = TableReader(os.path.join(simOutDir, "Main"))

		state_names = evaluationTime.readAttribute("state_names")
		process_names = evaluationTime.readAttribute("process_names")
		listener_names = evaluationTime.readAttribute("listener_names")
		logger_names = evaluationTime.readAttribute("logger_names")

		clock_times = evaluationTime.readColumn("clock_time")
		update_queries = evaluationTime.readColumn("update_queries_times")
		partition = evaluationTime.readColumn("partition_times")
		merge = evaluationTime.readColumn("merge_times")
		calculate_mass = evaluationTime.readColumn("calculate_mass_times")
		calculate_request = evaluationTime.readColumn("calculate_request_times")
		evolve_state = evaluationTime.readColumn("evolve_state_times")
		update = evaluationTime.readColumn("update_times")
		append = evaluationTime.readColumn("append_times")

		initialTime = mainReader.readAttribute("initialTime")
		time = (mainReader.readColumn("time") - initialTime) / 60  # min

		fig = plt.figure(figsize=(12, 15))
		gs = GridSpec(4, 2)
		subplot(gs[0, 0], time, update_queries, 'State.updateQueries', state_names)
		subplot(gs[1, 0], time, partition, 'State.partition', state_names)
		subplot(gs[2, 0], time, merge, 'State.merge', state_names)
		subplot(gs[3, 0], time, calculate_mass, 'State.calculateMass', state_names)
		subplot(gs[0, 1], time, calculate_request, 'Process.calculateRequest', process_names, sort=True)
		subplot(gs[1, 1], time, evolve_state, 'Process.evolveState', process_names, sort=True)
		subplot(gs[2, 1], time, update, 'Listener.update', listener_names, sort=True)
		subplot(gs[3, 1], time, append, 'Logger.append', logger_names)

		total_sim_eval_time = (clock_times[-1] - clock_times[0]) / 60  # min
		accounted_eval_time = np.sum(np.hstack((
			update_queries, partition, merge, calculate_mass,
			calculate_request, evolve_state, update, append))) / 60  # min

		fig.suptitle(
			'Total sim evaluation time = {:.2f} mins, Accounted evaluation time = {:.2f} mins'.format(
			total_sim_eval_time, accounted_eval_time))

		gs.tight_layout(fig)
		gs.update(top=0.95)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
