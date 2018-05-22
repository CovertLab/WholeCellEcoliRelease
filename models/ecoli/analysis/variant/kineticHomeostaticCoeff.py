
import argparse
import os
import re

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):
	if metadata["variant"] != "metabolismKineticHomeostaticRatio":
		print "This plot only runs for the 'metabolismKineticHomeostaticRatio' variant."
		return

	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"

	ap = AnalysisPaths(inputDir, variant_plot = True)
	variants = sorted(ap._path_data['variant'].tolist()) # Sorry for accessing private data

	if 0 in variants:
		variants.remove(0)

	if len(variants) == 0:
		return

	all_cells = sorted(ap.get_cells(variant = variants, seed = [0], generation = [0]))

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	barWidth = .35
	logMeanDifferences = []
	logStds = []
	for variant, simDir in zip(variants, all_cells):
		sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))
		simOutDir = os.path.join(simDir, "simOut")

		enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
		errorRelativeDifferences = enzymeKineticsdata.readColumn("kineticTargetRelativeDifferences")
		errorFluxNames = enzymeKineticsdata.readAttribute("kineticTargetFluxNames")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		enzymeKineticsdata.close()

		meanRelativeDiffTimeCourse = np.mean(np.abs(errorRelativeDifferences),axis=1)
		relativeDifferenceStd = np.std(meanRelativeDiffTimeCourse)

		logMeanDifferences.append(np.log10(meanRelativeDiffTimeCourse.mean()))
		logStds.append(np.log10(relativeDifferenceStd))


	plt.bar(np.arange(len(logMeanDifferences)), logMeanDifferences, barWidth, yerr=logStds)
	plt.xticks(np.arange(len(variants)), variants)
	plt.title("Kinetic/Homeostatic Ratio vs Deviation from Kinetic Target")
	plt.xlabel("Variant Number")
	plt.ylabel("Log10 Mean Absolute Relative Difference from Kinetic Targets")

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
