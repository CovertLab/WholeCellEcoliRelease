"""
Plots Figure 5D. Data is obtained from an analysis of compendium of environment shifts.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2017
"""

from __future__ import absolute_import

import os

import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

data = [1143, 330, 1920, 79, 854, 2, 25] #order is 100, 010, 110, 001, 101, 011, 111
dataArea = [1143, 330, 1920, 79, 854, 500, 25]


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		fig = plt.figure()
		v = venn3(subsets = dataArea, set_labels = ["0 < Freq. < 1", "Freq. = 1", "Freq. = 0"])
		ids = ["100", "010", "001", "110", "101", "011", "111"]

		for i, color in zip(ids, ["b", "r", "y", "purple", "green", "orange", "white"]):
			v.get_patch_by_id(i).set_color(color)

		for label, color in zip(v.set_labels, ["b", "r", "y"]):
			label.set_color(color)

		v.get_label_by_id("011").set_text("%s" % data[-2])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		for i in ids:
			v.get_label_by_id(i).set_text("")

		for label in v.set_labels:
			label.set_text("")

		plt.savefig(os.path.join(plotOutDir, plotOutFileName + "__clean.pdf"))
		plt.close()


if __name__ == "__main__":
	Plot().cli()
