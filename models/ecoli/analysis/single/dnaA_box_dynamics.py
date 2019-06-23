"""
Plot for dynamics of DnaA proteins binding to DnaA boxes across the chromosome

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/19
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os

from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		replication_data_reader = TableReader(os.path.join(simOutDir, 'ReplicationData'))

		# Load data
		initial_time = main_reader.readAttribute('initialTime')
		time = main_reader.readColumn('time') - initial_time
		free_DnaA_boxes = replication_data_reader.readColumn("free_DnaA_boxes")
		total_DnaA_boxes = replication_data_reader.readColumn("total_DnaA_boxes")
		criticalMassPerOriC = replication_data_reader.readColumn(
			"criticalMassPerOriC")

		names = [sim_data.moleculeIds.DnaA,
			sim_data.moleculeIds.DnaA_ATP_complex]
		(counts,) = read_bulk_molecule_counts(simOutDir, (names,))

		fig = plt.figure()
		fig.set_size_inches(8, 11.5)

		gs = gridspec.GridSpec(4, 1)

		# Plot critical cell mass per oriC
		ax = plt.subplot(gs[0, 0])
		ax.plot(time, criticalMassPerOriC)
		ax.plot(time, np.ones_like(time), "k--", linewidth=2)
		ax.set_yticks([0.5, 1.0])
		ax.set_ylabel("Critical mass\nper oriC")

		# Plot counts of free DnaA proteins and DnaA-ATP complexes
		ax = plt.subplot(gs[1, 0])
		ax.plot(time, counts[:, 0], label="Free DnaA proteins")
		ax.plot(time, counts[:, 1], label="Free DnaA-ATP complexes")
		ax.set_ylabel("Counts")
		ax.legend()

		# Plot counts of DnaA boxes
		ax = plt.subplot(gs[2, 0])
		ax.plot(time, total_DnaA_boxes, label="All DnaA boxes")
		ax.plot(time, free_DnaA_boxes, label="Free DnaA boxes")
		ax.set_ylabel("Counts")
		ax.legend()

		# Plot proportions of DnaA boxes bound to DnaA
		ax = plt.subplot(gs[3, 0])
		ax.plot(time, (total_DnaA_boxes - free_DnaA_boxes)/total_DnaA_boxes)
		ax.set_ylim([-0.05, 1.05])
		ax.set_xlabel("Time [s]")
		ax.set_ylabel("Proportions of DnaA boxes bound to DnaA")

		fig.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
