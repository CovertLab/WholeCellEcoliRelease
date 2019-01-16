"""
Produces histograms of frequency that production of a metabolite is limited (at least 50 time steps set by WINDOW)

@date: Created 1/12/2017
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

WINDOW = 50


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		sim_data = cPickle.load(open(simDataFile, "rb"))
		metaboliteNames = np.array(sorted(sim_data.process.metabolism.concDict.keys()))
		nMetabolites = len(metaboliteNames)

		fig, axesList = plt.subplots(3)
		fig.set_size_inches(11, 11)

		histo = np.zeros(4)
		limitedCounts = np.zeros(len(metaboliteNames))

		ax2 = axesList[2]
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			enzymeKineticsData = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
			metaboliteCounts = enzymeKineticsData.readColumn("metaboliteCountsFinal")
			normalizedCounts = metaboliteCounts / metaboliteCounts[1, :]
			enzymeKineticsData.close()

			# Read time info from the listener
			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			initialTime = main_reader.readAttribute("initialTime")
			time = main_reader.readColumn("time")

			metaboliteLimited = np.zeros((len(time), nMetabolites))

			diff = np.diff(normalizedCounts, axis = 0)
			limited = []
			for i in xrange(diff.shape[0] - WINDOW):
				currentStepLimited = np.where(np.any(diff[i:i + WINDOW] > 0, axis = 0) == False)[0].astype(int)
				metaboliteLimited[i, currentStepLimited] = 1
				limited = np.append(limited, currentStepLimited).astype(int)

			nLimited = len(np.unique(limited))
			if nLimited >= len(histo):
				histo = np.append(histo, np.zeros(nLimited - len(histo) + 1))
			histo[nLimited] += 1
			limitedCounts[limited] += 1

			ax2.plot(time / 60, metaboliteLimited * range(metaboliteLimited.shape[1]))
			ax2.axvline(initialTime / 60, color = "r", linestyle = "--")

		ax2.set_xlim([0, max(time) / 60])
		ax2.set_xlabel("Time (min)")
		ax2.set_ylabel("Limited")

		ax0 = axesList[0]
		labels = np.arange(len(histo))
		ax0.bar(labels - 0.5, histo, 1)
		ax0.set_xticks(labels)
		ax0.set_xlabel("Number of limited metabolites")
		ax0.set_ylabel("Number of generations")

		ax1 = axesList[1]
		ax1.bar(np.arange(len(np.where(limitedCounts > 0)[0])) - 0.4, limitedCounts[limitedCounts > 0])
		ax1.set_xticks(np.arange(len(np.where(limitedCounts > 0)[0])))
		ax1.set_xticklabels(metaboliteNames[limitedCounts > 0], fontsize = 6)
		ax1.set_xlabel("Metabolite Limited")
		ax1.set_ylabel("Number of genreations")

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")

def getMassData(simDir, massNames):
	simOutDir = os.path.join(simDir, "simOut")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
	mass = TableReader(os.path.join(simOutDir, "Mass"))

	massFractionData = np.zeros((len(massNames),time.size))

	for idx, massType in enumerate(massNames):
		massFractionData[idx,:] = mass.readColumn(massNames[idx])

	if len(massNames) == 1:
		massFractionData = massFractionData.reshape(-1)

	return time, massFractionData

def get_new_ylim(axis, new_min, new_max):
	ymin = axis.get_ylim()[0]
	ymax = axis.get_ylim()[1]

	if new_min < ymin:
		ymin = new_min
	if new_max > ymax:
		ymax = new_max

	return [ymin, ymax]

def removeNanReshape(a):
	return a[np.logical_not(np.isnan(a))]

if __name__ == "__main__":
	Plot().cli()
