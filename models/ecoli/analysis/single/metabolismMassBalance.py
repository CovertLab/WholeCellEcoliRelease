"""
Plot difference between mass imported to cell and mass created in metabolism process at each time step

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/27/2016
"""

from __future__ import absolute_import

import os

from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

		fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
		exFlux = fba_results.readColumn("externalExchangeFluxes")
		exMolec = fba_results.readAttribute("externalMoleculeIDs")
		fba_results.close()

		mass = TableReader(os.path.join(simOutDir, "Mass"))
		processMassDifferences = mass.readColumn("processMassDifferences")
		cellMass = mass.readColumn("dryMass")
		mass.close()

		exchangeMasses = {} # some duplicates in exMolec like CO2 and water so create dict to avoid double counting
		raw_data = KnowledgeBaseEcoli()
		for metabolite in raw_data.metabolites:
			for molecID in exMolec:
				if molecID.split("[")[0] == "WATER":
					exchangeMasses[molecID] = 18.015 * exFlux[:,exMolec.index(molecID)] * 10**-3 * cellMass * timeStepSec / 60 / 60
				if molecID.split("[")[0] == metabolite["id"]:
					exchangeMasses[molecID] = metabolite["mw7.2"] * exFlux[:,exMolec.index(molecID)] * 10**-3 * cellMass * timeStepSec / 60 / 60

		massInflux = 0
		for molecID in exchangeMasses.keys():
			massInflux += exchangeMasses[molecID]

		massProduced = processMassDifferences[:,0]	# in metabolism
		massDiff = massInflux + massProduced

		plt.plot(time / 60. / 60., -massDiff)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.ylabel("Mass Accumulation per time step (fg)")
		plt.title("Mass imported - mass created in metabolism process")

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
