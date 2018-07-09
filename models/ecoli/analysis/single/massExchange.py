"""
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/15/2016
"""

from __future__ import absolute_import
from __future__ import division

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FLUX_UNITS = "mmol/gDCW-hr"


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile))

		# Get exchange flux data
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		initialTime = units.s * TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = units.s * TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		externalExchangeFluxes = fbaResults.readColumn("externalExchangeFluxes")
		externalMoleculeIDs = np.array(fbaResults.readAttribute("externalMoleculeIDs"))
		fbaResults.close()

		massExchange = sim_data.getter.getMass(externalMoleculeIDs).asNumber(units.g / units.mmol) * externalExchangeFluxes # g / gDCW-hr

		# Get growth rate data
		growthRateData = TableReader(os.path.join(simOutDir, "Mass"))
		growthRate = ((1 / units.s) * growthRateData.readColumn("instantaniousGrowthRate")).asUnit(1 / units.h) # g / gDCW-hr
		doublingTime = (1 / growthRate) * np.log(2)

		# Plot stuff
		fig = plt.figure()
		fig.set_size_inches(8.5,11)

		ax1 = plt.subplot(3,1,1)
		ax1.plot(time.asNumber(units.min), doublingTime.asNumber(units.min))
		ax1.plot(time.asNumber(units.min), sim_data.doubling_time.asNumber(units.min) * np.ones(time.asNumber().size), linestyle='--')
		medianDoublingTime = np.median(doublingTime.asNumber(units.min)[1:])
		ax1.set_ylim([medianDoublingTime - 2*medianDoublingTime, medianDoublingTime + 2*medianDoublingTime])
		ax1.set_ylabel("Doubling\ntime (min)")

		ax2 = plt.subplot(3,1,2)
		ax2.plot(time.asNumber(units.min), massExchange)
		maxMassExchange = massExchange[100:].max()
		minMassExchange = massExchange[100:].min()
		ax2.set_ylim([minMassExchange, maxMassExchange])
		ax2.set_ylabel("Mass exchange\n(g / gDCW-hr)")

		ax3 = plt.subplot(3,1,3)
		water = massExchange[:, np.where(externalMoleculeIDs == "WATER[p]")[0][0]].copy()
		waterAll = massExchange[:, np.where(externalMoleculeIDs == "WATER[p]")[0][0]].copy()
		water[doublingTime.asNumber() > 0.] = np.nan
		ax3.plot(time.asNumber(units.min), waterAll, 'k.')
		ax3.plot(time.asNumber(units.min), water, 'b.')
		maxMassExchange = massExchange[100:].max()
		minMassExchange = massExchange[100:].min()
		ax3.set_ylim([minMassExchange, maxMassExchange])
		ax3.set_ylabel("Water exchange\nwhen doubling time < 0")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
