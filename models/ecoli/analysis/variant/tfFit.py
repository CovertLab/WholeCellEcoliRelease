from __future__ import absolute_import


import os

import numpy as np
from matplotlib import pyplot as plt
import bokeh.io
from bokeh.io import vplot
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import (HoverTool, BoxZoomTool, LassoSelectTool, PanTool,
	WheelZoomTool, ResizeTool, UndoTool, RedoTool)
from bokeh.models import CustomJS
from bokeh.models.widgets import Button

import cPickle
import scipy.stats

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

NUMERICAL_ZERO = 1e-12


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata["variant"] != "tf_activity":
			print "This plot only runs for the 'tf_activity' variant."
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
		targetId = []
		targetCondition = []
		targetToTfType = {}

		for variant, simDir in zip(variants, all_cells):
			sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))

			delta_prob = sim_data.process.transcription_regulation.delta_prob

			tfList = ["basal (no TF)"] + sorted(sim_data.tfToActiveInactiveConds)
			simOutDir = os.path.join(simDir, "simOut")
			tf = tfList[(variant + 1) // 2]

			if variant % 2 == 1:
				tfStatus = "active"
			else:
				tfStatus = "inactive"

			rna_synth_prob_reader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
			rna_ids = rna_synth_prob_reader.readAttribute("rnaIds")
			tf_ids = rna_synth_prob_reader.readAttribute("tf_ids")
			n_bound_TF_per_TU = rna_synth_prob_reader.readColumn(
				"n_bound_TF_per_TU").reshape((-1, len(rna_ids), len(tf_ids)))
			gene_copy_number = rna_synth_prob_reader.readColumn("gene_copy_number")

			tf_idx = tf_ids.index(tf)
			tf_targets = sim_data.tfToFC[tf]
			tf_target_indexes = np.array([
				rna_ids.index(tf_target + "[c]") for tf_target in tf_targets])

			tfTargetBoundCountsAll = n_bound_TF_per_TU[:, tf_target_indexes, tf_idx]
			tfTargetSynthProbAll = rna_synth_prob_reader.readColumn("rnaSynthProb")[:, tf_target_indexes]
			tf_target_gene_copies_all = gene_copy_number[:, tf_target_indexes]

			for i, tfTarget in enumerate(sorted(sim_data.tfToFC[tf])):
				tfTargetBoundCounts = tfTargetBoundCountsAll[:, i].reshape(-1)
				tf_target_copies = tf_target_gene_copies_all[:, i].reshape(-1)

				expectedProbBound.append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
				simulatedProbBound.append(
					(tfTargetBoundCounts[5:].astype(np.float64)/tf_target_copies[5:]).mean())

				tfTargetSynthProb = tfTargetSynthProbAll[:, i].reshape(-1)

				target_idx = rna_ids.index(tfTarget + "[c]")
				regulating_tf_idxs = delta_prob['deltaJ'][delta_prob['deltaI'] == target_idx]

				for i in regulating_tf_idxs:
					if tfTarget not in targetToTfType:
						targetToTfType[tfTarget] = []
					targetToTfType[tfTarget].append(sim_data.process.transcription_regulation.tfToTfType[tf_ids[i]])

				expectedSynthProb.append(sim_data.process.transcription.rnaSynthProb[tf + "__" + tfStatus][target_idx])
				simulatedSynthProb.append(tfTargetSynthProb[5:].mean())

				targetId.append(tfTarget)
				targetCondition.append(tf + "__" + tfStatus)


		expectedProbBound = np.array(expectedProbBound)
		simulatedProbBound = np.array(simulatedProbBound)
		expectedSynthProb = np.array(expectedSynthProb)
		simulatedSynthProb = np.array(simulatedSynthProb)

		regressionResult = scipy.stats.linregress(np.log10(expectedProbBound[expectedProbBound > NUMERICAL_ZERO]), np.log10(simulatedProbBound[expectedProbBound > NUMERICAL_ZERO]))
		regressionResultLargeValues = scipy.stats.linregress(np.log10(expectedProbBound[expectedProbBound > 1e-2]), np.log10(simulatedProbBound[expectedProbBound > 1e-2]))

		ax = plt.subplot(2, 1, 1)
		ax.scatter(np.log10(expectedProbBound), np.log10(simulatedProbBound))
		plt.xlabel("log10(Expected probability bound)", fontsize = 6)
		plt.ylabel("log10(Simulated probability bound)", fontsize = 6)
		plt.title("Slope: %0.3f   Intercept: %0.3e      (Without Small Values:  Slope: %0.3f Intercept: %0.3e)" % (regressionResult.slope, regressionResult.intercept, regressionResultLargeValues.slope, regressionResultLargeValues.intercept), fontsize = 6)
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

		regressionResult = scipy.stats.linregress(np.log10(expectedSynthProb[expectedSynthProb > NUMERICAL_ZERO]), np.log10(simulatedSynthProb[expectedSynthProb > NUMERICAL_ZERO]))

		ax = plt.subplot(2, 1, 2)
		ax.scatter(np.log10(expectedSynthProb), np.log10(simulatedSynthProb))
		plt.xlabel("log10(Expected synthesis probability)", fontsize = 6)
		plt.ylabel("log10(Simulated synthesis probability)", fontsize = 6)
		plt.title("Slope: %0.3f   Intercept: %0.3e" % (regressionResult.slope, regressionResult.intercept), fontsize = 6)
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		# Probability bound - hover for ID
		source1 = ColumnDataSource(data = dict(x = np.log10(expectedProbBound), y = np.log10(simulatedProbBound), ID = targetId, condition = targetCondition))
		hover1 = HoverTool(tooltips = [("ID", "@ID"), ("condition", "@condition")])
		tools1 = [hover1, BoxZoomTool(), LassoSelectTool(), PanTool(), WheelZoomTool(), ResizeTool(),	UndoTool(),	RedoTool(), "reset"]
		s1 = figure(
			x_axis_label = "log10(Expected probability bound)",
			y_axis_label = "log10(Simulated probability bound)",
			width = 800,
			height = 500,
			tools = tools1)
		s1.scatter("x", "y", source = source1)

		if not os.path.exists(os.path.join(plotOutDir, "html_plots")):
			os.makedirs(os.path.join(plotOutDir, "html_plots"))
		bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + "__probBound" + ".html"), title = plotOutFileName, autosave = False)
		bokeh.io.save(s1)

		# Synthesis probability - hover for ID
		source2 = ColumnDataSource(data = dict(x = np.log10(expectedSynthProb), y = np.log10(simulatedSynthProb), ID = targetId, condition = targetCondition))
		hover2 = HoverTool(tooltips = [("ID", "@ID"), ("condition", "@condition")])
		tools2 = [hover2, BoxZoomTool(), LassoSelectTool(), PanTool(), WheelZoomTool(), ResizeTool(),	UndoTool(),	RedoTool(), "reset"]
		s2 = figure(
			x_axis_label = "log10(Expected synthesis probability)",
			y_axis_label = "log10(Simulated synthesis probability)",
			width = 800,
			height = 500,
			tools = tools2)
		s2.scatter("x", "y", source = source2)

		bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + "__synthProb" + ".html"), title = plotOutFileName, autosave = False)
		bokeh.io.save(s2)

		# Synthesis probability - filter targets by TF type
		bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + "__synthProb__interactive" + ".html"), title = plotOutFileName, autosave = False)

		tfTypes = []
		for i in targetId:
			if i in targetToTfType:
				uniqueSet = np.unique(targetToTfType[i])

				if uniqueSet.shape[0] == 1:
					tfTypes.append(uniqueSet[0])
				elif uniqueSet.shape[0] == 3:
					tfTypes.append("all")
				else:
					tfTypes.append(uniqueSet[0] + "_" + uniqueSet[1])
			else:
				tfTypes.append("none")
		tfTypes = np.array(tfTypes)


		x0 = np.copy(expectedSynthProb)
		x0[np.where(tfTypes != "0CS")]= np.nan
		x1 = np.copy(expectedSynthProb)
		x1[np.where(tfTypes != "1CS")]= np.nan
		x2 = np.copy(expectedSynthProb)
		x2[np.where(tfTypes != "2CS")]= np.nan
		x01 = np.copy(expectedSynthProb)
		x01[np.where(tfTypes != "0CS_1CS")] = np.nan
		x02 = np.copy(expectedSynthProb)
		x02[np.where(tfTypes != "0CS_2CS")] = np.nan
		x12 = np.copy(expectedSynthProb)
		x12[np.where(tfTypes != "1CS_2CS")] = np.nan

		y0 = np.copy(simulatedSynthProb)
		y0[np.where(tfTypes != "0CS")] = np.nan
		y1 = np.copy(simulatedSynthProb)
		y1[np.where(tfTypes != "1CS")] = np.nan
		y2 = np.copy(simulatedSynthProb)
		y2[np.where(tfTypes != "2CS")] = np.nan
		y01 = np.copy(simulatedSynthProb)
		y01[np.where(tfTypes != "0CS_1CS")] = np.nan
		y02 = np.copy(simulatedSynthProb)
		y02[np.where(tfTypes != "0CS_2CS")] = np.nan
		y12 = np.copy(simulatedSynthProb)
		x12[np.where(tfTypes != "1CS_2CS")] = np.nan

		source_all = ColumnDataSource(data = dict(x = np.log10(expectedSynthProb), y = np.log10(simulatedSynthProb), ID = targetId, condition = targetCondition))
		source_tf = ColumnDataSource(data = dict(x0 = np.log10(x0), y0 = np.log10(y0), x1 = np.log10(x1), y1 = np.log10(y1), x2 = np.log10(x2), y2 = np.log10(y2), x01 = np.log10(x01), y01 = np.log10(y01), x02 = np.log10(x02), y02 = np.log10(y02), x12 = np.log10(x12), y12 = np.log10(y12),x123 = np.log10(expectedSynthProb), y123 = np.log10(simulatedSynthProb),ID = targetId, condition = targetCondition))
		hover3 = HoverTool(tooltips = [("ID", "@ID"), ("condition", "@condition")])
		tools3 = [hover3, BoxZoomTool(), LassoSelectTool(), PanTool(), WheelZoomTool(), ResizeTool(),	UndoTool(),	RedoTool(), "reset"]

		axis_max = np.ceil(np.log10(expectedSynthProb).max())
		for i in np.sort(expectedSynthProb):
			if i > 0:
				break
		axis_min = np.floor(np.log10(i))
		s3 = figure(
			x_axis_label = "log10(Expected synthesis probability)",
			y_axis_label = "log10(Simulated synthesis probability)",
			plot_width=800, plot_height=500,
			x_range = (axis_min, axis_max),
			y_range = (axis_min, axis_max),
			tools = tools3,
			)
		s3.scatter("x", "y", source = source_all)
		callback = CustomJS(args = dict(source_all = source_all, source_tf = source_tf), code =
			"""
			var data_all = source_all.get('data');
			var data_tf = source_tf.get('data');
			data_all['x'] = data_tf['x' + cb_obj.get("name")];
			data_all['y'] = data_tf['y' + cb_obj.get("name")];
			source_all.trigger('change');
			""")

		toggle0 = Button(label = "0CS", callback = callback, name = "0")
		toggle1 = Button(label = "1CS", callback = callback, name = "1")
		toggle2 = Button(label = "2CS", callback = callback, name = "2")
		toggle3 = Button(label = "0CS and 1CS", callback = callback, name = "01")
		toggle4 = Button(label = "0CS and 2CS", callback = callback, name = "02")
		toggle5 = Button(label = "1CS and 2CS", callback = callback, name = "12")
		toggle6 = Button(label = "All", callback = callback, name = "123")
		layout = vplot(toggle0, toggle1, toggle2, toggle3, toggle4, toggle5, toggle6, s3)
		bokeh.io.save(layout)
		bokeh.io.curstate().reset()


if __name__ == "__main__":
	Plot().cli()
