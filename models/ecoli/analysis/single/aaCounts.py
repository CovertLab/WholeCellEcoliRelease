"""
Plot amino acid counts

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import absolute_import, division, print_function

import os
import cPickle

from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile))

		aaIDs = sim_data.moleculeGroups.aaIDs
		(aaCounts,) = read_bulk_molecule_counts(simOutDir, (aaIDs,))

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		plt.figure(figsize = (8.5, 11))

		for idx in xrange(21):

			plt.subplot(6, 4, idx + 1)

			plt.plot(time / 60., aaCounts[:, idx], linewidth = 2)
			plt.xlabel("Time (min)")
			plt.ylabel("Counts")
			plt.title(aaIDs[idx], fontsize=8)
			plt.tick_params(labelsize=8)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
