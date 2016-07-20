
import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile, metadata = None):
	if metadata["variant"] != "tfActivity":
		print "This plot only runs for the 'tfActivity' variant."
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

	expectedProbBound = []
	simulatedProbBound = []
	expectedSynthProb = []
	simulatedSynthProb = []


	for variant, simDir in zip(variants, all_cells):
		sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))
		tfList = ["basal (no TF)"] + sorted(sim_data.tfToActiveInactiveConds)
		simOutDir = os.path.join(simDir, "simOut")
		tf = tfList[(variant + 1) // 2]
		tfStatus = None
		if variant % 2 == 1:
			tfStatus = "active"
		else:
			tfStatus = "inactive"

		bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")

		rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
		rnaIds = rnaSynthProbReader.readAttribute("rnaIds")


		for tfTarget in sorted(sim_data.tfToFC[tf]):

			tfTargetBoundId = [tfTarget + "__" + tf]
			tfTargetBoundIndex = np.array([bulkMoleculeIds.index(x) for x in tfTargetBoundId])
			tfTargetBoundCounts = bulkMoleculesReader.readColumn("counts")[:, tfTargetBoundIndex].reshape(-1)

			expectedProbBound.append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
			simulatedProbBound.append(tfTargetBoundCounts.mean())


			tfTargetSynthProbId = [tfTarget + "[c]"]
			tfTargetSynthProbIndex = np.array([rnaIds.index(x) for x in tfTargetSynthProbId])
			tfTargetSynthProb = rnaSynthProbReader.readColumn("rnaSynthProb")[:, tfTargetSynthProbIndex].reshape(-1)

			rnaIdx = np.where(sim_data.process.transcription.rnaData["id"] == tfTarget + "[c]")[0][0]

			expectedSynthProb.append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][rnaIdx])
			simulatedSynthProb.append(tfTargetSynthProb.mean())

		bulkMoleculesReader.close()
		rnaSynthProbReader.close()


	ax = plt.subplot(2, 1, 1)
	ax.scatter(expectedProbBound, simulatedProbBound)
	plt.xlabel("Expected probability bound", fontsize = 6)
	plt.ylabel("Simulated probability bound", fontsize = 6)
	ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

	ax = plt.subplot(2, 1, 2)
	ax.scatter(np.log10(expectedSynthProb), np.log10(simulatedSynthProb))
	plt.xlabel("log10(Expected synthesis probability)", fontsize = 6)
	plt.ylabel("log10(Simulated synthesis probability)", fontsize = 6)
	ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

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
