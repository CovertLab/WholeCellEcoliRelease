"""
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/8/2014
"""

from __future__ import absolute_import, division, print_function

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

GLUCOSE_ID = "GLC[p]"

FLUX_UNITS = units.mmol / units.g / units.h
MASS_UNITS = units.fg
GROWTH_UNITS = units.fg / units.s


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile, "rb"))

		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		externalExchangeFluxes = fbaResults.readColumn("externalExchangeFluxes")

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime
		timeStepSec = main_reader.readColumn("timeStepSec")

		externalMoleculeIDs = np.array(fbaResults.readAttribute("externalMoleculeIDs"))
		fbaResults.close()

		if GLUCOSE_ID not in externalMoleculeIDs:
			print("This plot only runs when glucose is the carbon source.")
			return

		glucoseIdx = np.where(externalMoleculeIDs == GLUCOSE_ID)[0][0]
		glucoseFlux = FLUX_UNITS * externalExchangeFluxes[:, glucoseIdx]

		mass = TableReader(os.path.join(simOutDir, "Mass"))
		# cellMass = MASS_UNITS * mass.readColumn("cellMass")
		cellDryMass = MASS_UNITS * mass.readColumn("dryMass")
		growth = GROWTH_UNITS * mass.readColumn("growth") / timeStepSec
		mass.close()

		glucoseMW = sim_data.getter.getMass([GLUCOSE_ID])[0]

		glucoseMassFlux = glucoseFlux * glucoseMW * cellDryMass

		glucoseMassYield = growth / -glucoseMassFlux

		fig = plt.figure(figsize = (8.5, 11))
		plt.plot(time, glucoseMassYield)
		plt.xlabel("Time (s)")
		plt.ylabel("g cell / g glucose")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
