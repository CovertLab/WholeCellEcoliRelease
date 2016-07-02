#!/usr/bin/env python
"""
Plot difference between mass imported to cell and mass created in metabolism process at each time step

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/27/2016
"""

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli



def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata=None):
	
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

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
