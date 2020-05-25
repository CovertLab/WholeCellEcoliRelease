"""
Plot trp regulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/2016
"""

from __future__ import absolute_import, division, print_function

import os

from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

BURN_IN = 10


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		allDirs = ap.get_cells()

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		trpIdx = sim_data.moleculeGroups.aaIDs.index("TRP[c]")

		plt.figure(figsize = (8.5, 11))

		for simDir in allDirs:
			simOutDir = os.path.join(simDir, "simOut")

			growthLimits = TableReader(os.path.join(simOutDir, "GrowthLimits"))

			trpRequests = growthLimits.readColumn("aaRequestSize")[BURN_IN:, trpIdx]

			growthLimits.close()

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

			moleculeIds = bulkMolecules.readAttribute("objectNames")

			trpSynIdx = moleculeIds.index("TRYPSYN[c]")

			trpSynCounts = bulkMolecules.readColumn("counts")[BURN_IN:, trpSynIdx]

			bulkMolecules.close()

			trpSynKcat = 2**( (37. - 25.) / 10.) * 4.1 # From PMID 6402362 (kcat of 4.1/s measured at 25 C)

			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			# initialTime = main_reader.readAttribute("initialTime")
			time = main_reader.readColumn("time")[BURN_IN:]
			timeStep = main_reader.readColumn("timeStepSec")[BURN_IN:]


			trpSynMaxCapacity = trpSynKcat * trpSynCounts * timeStep


			##############################################################
			ax = self.subplot(3, 1, 1)
			ax.plot(time / 60., trpSynMaxCapacity, color = "b")
			plt.ylabel("Tryptophan Synthase Max Capacity", fontsize = 10)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 10)
			ax.set_xticks([])
			##############################################################

			##############################################################
			ax = self.subplot(3, 1, 2)
			ax.plot(time, trpRequests, color = "b")
			plt.ylabel("Trp Requested By Translation", fontsize = 10)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize=10)
			ax.set_xticks([])
			##############################################################


			##############################################################
			ax = self.subplot(3, 1, 3)
			ax.plot(time / 3600., trpSynMaxCapacity / trpRequests, color = "b")
			ax.plot([0, time[-1] / 3600.], [1., 1.], "k--")
			plt.ylabel("(Max capacity) / (Request)", fontsize = 10)

			ymin, ymax = ax.get_ylim()
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2f" % ymin, "%0.2f" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			# ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 10)
			ax.set_xticks(ax.get_xlim())
			##############################################################

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
