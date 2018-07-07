from __future__ import absolute_import


import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle
import scipy.stats

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

NUMERICAL_ZERO = 1e-12


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
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

		expectedProbBound = [[], [], []]
		simulatedProbBound = [[], [], []]
		expectedSynthProb = [[], [], []]
		simulatedSynthProb = [[], [], []]

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

			tfTargetBoundIds = []
			tfTargetBoundIndices = []
			tfTargetSynthProbIds = []
			tfTargetSynthProbIndices = []
			for tfTarget in sorted(sim_data.tfToFC[tf]):
				tfTargetBoundIds.append(tfTarget + "__" + tf)
				tfTargetBoundIndices.append(bulkMoleculeIds.index(tfTargetBoundIds[-1]))
				tfTargetSynthProbIds.append(tfTarget + "[c]")
				tfTargetSynthProbIndices.append(rnaIds.index(tfTargetSynthProbIds[-1]))
			tfTargetBoundCountsAll = bulkMoleculesReader.readColumn("counts")[:, tfTargetBoundIndices]
			tfTargetSynthProbAll = rnaSynthProbReader.readColumn("rnaSynthProb")[:, tfTargetSynthProbIndices]

			for targetIdx, tfTarget in enumerate(sorted(sim_data.tfToFC[tf])):
				tfTargetBoundCounts = tfTargetBoundCountsAll[:, targetIdx].reshape(-1)

				tfTargetSynthProbId = [tfTarget + "[c]"]
				tfTargetSynthProbIndex = np.array([rnaIds.index(x) for x in tfTargetSynthProbId])
				tfTargetSynthProb = rnaSynthProbReader.readColumn("rnaSynthProb")[:, tfTargetSynthProbIndex].reshape(-1)

				rnaIdx = np.where(sim_data.process.transcription.rnaData["id"] == tfTarget + "[c]")[0][0]

				tfType = sim_data.process.transcription_regulation.tfToTfType[tf]
				if tfType == "0CS":
					expectedProbBound[0].append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
					simulatedProbBound[0].append(tfTargetBoundCounts[5:].mean())

					expectedSynthProb[0].append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][rnaIdx])
					simulatedSynthProb[0].append(tfTargetSynthProb[5:].mean())

				elif tfType == "1CS":
					expectedProbBound[1].append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
					simulatedProbBound[1].append(tfTargetBoundCounts[5:].mean())

					expectedSynthProb[1].append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][rnaIdx])
					simulatedSynthProb[1].append(tfTargetSynthProb[5:].mean())
				else:
					expectedProbBound[2].append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
					simulatedProbBound[2].append(tfTargetBoundCounts[5:].mean())

					expectedSynthProb[2].append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][rnaIdx])
					simulatedSynthProb[2].append(tfTargetSynthProb[5:].mean())

			bulkMoleculesReader.close()
			rnaSynthProbReader.close()

		expectedProbBound = np.array(expectedProbBound)
		simulatedProbBound = np.array(simulatedProbBound)
		expectedSynthProb = np.array(expectedSynthProb)
		simulatedSynthProb = np.array(simulatedSynthProb)

		rows = 2
		cols = 3

		for i in np.arange(3):
			iExpectedProbBound = np.array(expectedProbBound[i])
			iSimulatedProbBound = np.array(simulatedProbBound[i])
			iExpectedSynthProb = np.array(expectedSynthProb[i])
			iSimulatedSynthProb = np.array(simulatedSynthProb[i])

			regressionResult = scipy.stats.linregress(np.log10(iExpectedProbBound[iExpectedProbBound > NUMERICAL_ZERO]), np.log10(iSimulatedProbBound[iExpectedProbBound > NUMERICAL_ZERO]))
			regressionResultLargeValues = scipy.stats.linregress(np.log10(iExpectedProbBound[iExpectedProbBound > 1e-2]), np.log10(iSimulatedProbBound[iExpectedProbBound > 1e-2]))

			ax = plt.subplot(rows, cols, i + 1)
			ax.scatter(np.log10(expectedProbBound[i]), np.log10(simulatedProbBound[i]))
			ax.set_xlim([-3.0, 0.5])
			ax.set_ylim([-2.5, 0.5])
			plt.xlabel("log10(Expected probability bound)", fontsize = 6)
			plt.ylabel("log10(Simulated probability bound)", fontsize = 6)
			plt.title("%sCS\nSlope: %0.3f   Intercept: %0.3e      \n(Without Small Values:\nSlope: %0.3f Intercept: %0.3e)" % (i, regressionResult.slope, regressionResult.intercept, regressionResultLargeValues.slope, regressionResultLargeValues.intercept), fontsize = 6)
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

			regressionResult = scipy.stats.linregress(np.log10(iExpectedSynthProb[iExpectedSynthProb > NUMERICAL_ZERO]), np.log10(iSimulatedSynthProb[iExpectedSynthProb > NUMERICAL_ZERO]))

			ax = plt.subplot(rows, cols, i + 4)
			ax.scatter(np.log10(expectedSynthProb[i]), np.log10(simulatedSynthProb[i]))
			plt.xlabel("log10(Expected synthesis probability)", fontsize = 6)
			plt.ylabel("log10(Simulated synthesis probability)", fontsize = 6)
			plt.title("Slope: %0.3f   Intercept: %0.3e" % (regressionResult.slope, regressionResult.intercept), fontsize = 6)
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

		plt.subplots_adjust(hspace = 0.4, wspace = 0.4)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		rows = 2
		cols = 1
		color = ["blue", "green", "red"]

		for i in np.arange(3):
			iExpectedProbBound = np.array(expectedProbBound[i])
			iSimulatedProbBound = np.array(simulatedProbBound[i])
			iExpectedSynthProb = np.array(expectedSynthProb[i])
			iSimulatedSynthProb = np.array(simulatedSynthProb[i])

			ax = plt.subplot(rows, cols, 1)
			ax.scatter(np.log10(expectedProbBound[i]), np.log10(simulatedProbBound[i]), color = color[i], alpha = 0.5)
			plt.xlabel("log10(Expected probability bound)", fontsize = 6)
			plt.ylabel("log10(Simulated probability bound)", fontsize = 6)
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

			ax = plt.subplot(rows, cols, 2)
			ax.scatter(np.log10(expectedSynthProb[i]), np.log10(simulatedSynthProb[i]), color = color[i], alpha = 0.5)
			plt.xlabel("log10(Expected synthesis probability)", fontsize = 6)
			plt.ylabel("log10(Simulated synthesis probability)", fontsize = 6)
			ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

		plt.subplots_adjust(hspace = 0.4, wspace = 0.4)
		exportFigure(plt, plotOutDir, plotOutFileName + "__overlap", metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
