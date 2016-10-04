#!/usr/bin/env python
"""
Plot objective reaction components over time

@date: Created 9/27/2016
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

NUMERICAL_ZERO = 1e-20

CLOSENESS_THRESHOLD = 1

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaData = TableReader(os.path.join(simOutDir, "FBAResults"))
	
	homeostaticTargetMolecules = fbaData.readAttribute("homeostaticTargetMolecules")
	kineticTargetFluxNames = fbaData.readAttribute("kineticTargetFluxNames")

	homeostaticObjectiveComponents = fbaData.readColumn("homeostaticObjectiveValues")
	kineticObjectiveComponents = fbaData.readColumn("kineticObjectiveValues")
	homeostaticObjectiveWeights = fbaData.readColumn("homeostaticObjectiveWeight")
	
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	fbaData.close()

	plt.figure(figsize = (8.5, 11))
	plt.suptitle("Metabolism Objective Components")

	plt.subplot(3,1,1)

	plt.plot(time / 60, homeostaticObjectiveWeights * np.sum(homeostaticObjectiveComponents,axis=1), label="Homeostatic Objective Component")

	plt.plot(time / 60, (1-homeostaticObjectiveWeights) * np.sum(kineticObjectiveComponents,axis=1), label="Kinetic Objective Component")

	plt.xlabel("Time (min)")
	plt.ylabel("Sum of Objective components", fontsize="x-small")
	plt.legend(framealpha=.5, fontsize=6)


	plt.subplot(3,1,2)

	plt.plot(time / 60, np.mean(homeostaticObjectiveComponents < NUMERICAL_ZERO, axis=1), label="Homeostatic Objective Components")

	plt.plot(time / 60, np.mean(kineticObjectiveComponents < NUMERICAL_ZERO, axis=1), label="Kinetic Objective Components")

	plt.xlabel("Time (min)")
	plt.ylabel("Portion of objective components less than numerical zero" ,fontsize="x-small")
	plt.legend(framealpha=.5, fontsize=6)

	plt.subplot(3,1,3)

	plt.plot(time / 60, np.mean(homeostaticObjectiveComponents < CLOSENESS_THRESHOLD, axis=1), label="Homeostatic Objective Components")

	plt.plot(time / 60, np.mean(kineticObjectiveComponents < CLOSENESS_THRESHOLD, axis=1), label="Kinetic Objective Components")

	plt.xlabel("Time (min)")
	plt.ylabel("Portion of objective components less than {}".format(CLOSENESS_THRESHOLD) ,fontsize="x-small")
	plt.legend(framealpha=.5, fontsize=6)


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
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
