from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

NUMERICAL_ZERO = 1e-12


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata["variant"] != "tf_activity":
			print("This plot only runs for the 'tf_activity' variant.")
			return

		ap = AnalysisPaths(inputDir, variant_plot = True)
		variants = sorted(ap._path_data['variant'].tolist()) # Sorry for accessing private data

		if 0 in variants:
			variants.remove(0)

		if len(variants) == 0:
			return

		all_cells = sorted(ap.get_cells(variant = variants, seed = [0], generation = [0]))

		expectedProbBound = [[], [], []]
		simulatedProbBound = [[], [], []]
		expectedSynthProb = [[], [], []]
		simulatedSynthProb = [[], [], []]

		for variant, simDir in zip(variants, all_cells):

			sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))
			tfList = ["basal (no TF)"] + sorted(sim_data.tfToActiveInactiveConds)
			simOutDir = os.path.join(simDir, "simOut")
			tf = tfList[(variant + 1) // 2]

			if variant % 2 == 1:
				tfStatus = "active"
			else:
				tfStatus = "inactive"

			rna_synth_prob_reader = TableReader(
				os.path.join(simOutDir, "RnaSynthProb"))
			rna_ids = rna_synth_prob_reader.readAttribute("rnaIds")
			tf_ids = rna_synth_prob_reader.readAttribute("tf_ids")
			n_bound_TF_per_TU = rna_synth_prob_reader.readColumn(
				"n_bound_TF_per_TU").reshape(
				(-1, len(rna_ids), len(tf_ids)))
			gene_copy_number = rna_synth_prob_reader.readColumn(
				"gene_copy_number")

			tf_idx = tf_ids.index(tf)
			tf_targets = sim_data.tfToFC[tf]
			tf_target_indexes = np.array([
				rna_ids.index(tf_target + "[c]") for tf_target in tf_targets
				])

			tfTargetBoundCountsAll = n_bound_TF_per_TU[:, tf_target_indexes, tf_idx]
			tfTargetSynthProbAll = rna_synth_prob_reader.readColumn("rnaSynthProb")[:, tf_target_indexes]
			tf_target_gene_copies_all = gene_copy_number[:, tf_target_indexes]

			for i, tfTarget in enumerate(sorted(sim_data.tfToFC[tf])):
				tfTargetBoundCounts = tfTargetBoundCountsAll[:, i].reshape(-1)
				tf_target_copies = tf_target_gene_copies_all[:, i].reshape(-1)
				tfTargetSynthProb = tfTargetSynthProbAll[:, i].reshape(-1)

				rnaIdx = np.where(sim_data.process.transcription.rnaData["id"] == tfTarget + "[c]")[0][0]

				tfType = sim_data.process.transcription_regulation.tfToTfType[tf]

				if tfType == "0CS":
					expectedProbBound[0].append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
					simulatedProbBound[0].append(
						(tfTargetBoundCounts[5:].astype(np.float64)/tf_target_copies[5:]).mean())

					expectedSynthProb[0].append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][rnaIdx])
					simulatedSynthProb[0].append(tfTargetSynthProb[5:].mean())

				elif tfType == "1CS":
					expectedProbBound[1].append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
					simulatedProbBound[1].append(
						(tfTargetBoundCounts[5:].astype(np.float64)/tf_target_copies[5:]).mean())

					expectedSynthProb[1].append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][rnaIdx])
					simulatedSynthProb[1].append(tfTargetSynthProb[5:].mean())
				else:
					expectedProbBound[2].append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
					simulatedProbBound[2].append(
						(tfTargetBoundCounts[5:].astype(np.float64)/tf_target_copies[5:]).mean())

					expectedSynthProb[2].append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][rnaIdx])
					simulatedSynthProb[2].append(tfTargetSynthProb[5:].mean())

		expectedProbBound = np.array(expectedProbBound)
		simulatedProbBound = np.array(simulatedProbBound)
		expectedSynthProb = np.array(expectedSynthProb)
		simulatedSynthProb = np.array(simulatedSynthProb)

		rows = 2
		cols = 1
		color = ["blue", "green", "red"]

		for i in np.arange(3):
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
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
