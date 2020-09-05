"""
Plot NTP counts

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile, 'rb'))

		dntpIDs = sim_data.molecule_groups.dntps
		(dntpCounts,) = read_bulk_molecule_counts(simOutDir, (dntpIDs,))

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		plt.figure(figsize = (8.5, 8.5))

		for idx in range(4):

			plt.subplot(2, 2, idx + 1)

			plt.plot(time / 60., dntpCounts[:, idx], linewidth = 2)
			plt.xlabel("Time (min)")
			plt.ylabel("Counts")
			plt.title(dntpIDs[idx])

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
