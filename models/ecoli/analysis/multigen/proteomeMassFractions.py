"""
Plot normalized mass fractions of proteome for each protein monomer

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/12/16
"""

from __future__ import absolute_import, division

import os
import cPickle
import sys

import numpy as np
from matplotlib import pyplot as plt
import scipy.cluster

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		allDirs = ap.get_cells()

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		ids_translation = sim_data.process.translation.monomerData["id"].tolist()

		proteomeMWs = sim_data.getter.getMass(ids_translation)

		time = []
		proteomeMassFractions = []
		initialTime = 0

		plt.figure(figsize = (8.5, 11))

		for simDir in allDirs:
			simOutDir = os.path.join(simDir, "simOut")

			# Load time
			time = np.append(time, TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime)

			# Load mass data
			massReader = TableReader(os.path.join(simOutDir, "Mass"))
			proteinMass = units.fg * massReader.readColumn("proteinMass")
			massReader.close()

			monomerCounts = TableReader(os.path.join(simOutDir, "MonomerCounts"))
			counts = monomerCounts.readColumn("monomerCounts")

			# Get mass of proteins in cell
			proteomeMasses = 1. / nAvogadro * counts * proteomeMWs
			if len(proteomeMassFractions):
				proteomeMassFractions = np.concatenate((proteomeMassFractions, proteomeMasses.asNumber(units.fg).T / proteinMass.asNumber(units.fg)), axis = 1)
			else:
				proteomeMassFractions = proteomeMasses.asNumber(units.fg).T / proteinMass.asNumber(units.fg)

		# Prevent divide by 0
		proteomeMassFractions += 1e-9

		# Normalize based on mass at early time point and cluster similar traces
		normalizedMassFractions = proteomeMassFractions.T / proteomeMassFractions.T[10, :]
		sys.setrecursionlimit(10000)
		linkage = scipy.cluster.hierarchy.linkage(normalizedMassFractions.T)
		dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot = True)

		nPlots = 16
		nProteins = normalizedMassFractions.shape[1]

		for i in range(nPlots):
			ax = plt.subplot(nPlots, 1, i + 1)
			indStart = i * nProteins // nPlots
			indEnd = (i + 1) * nProteins // nPlots
			indexes = dendro["leaves"][indStart:indEnd]

			proteomeSubset = normalizedMassFractions[:, indexes]

			ax.plot(time, proteomeSubset)
			plt.ylabel("Mass Fraction of Proteome", fontsize = 6)

			ymin = np.amin(proteomeSubset * 0.9)
			ymax = np.amax(proteomeSubset * 1.1)
			if ymin != ymax:
				ax.set_ylim([ymin, ymax])
				ax.set_yticks([ymin, ymax])
				ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
			ax.set_xticks([])

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
