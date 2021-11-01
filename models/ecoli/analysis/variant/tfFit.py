from __future__ import absolute_import, division, print_function

import os

import bokeh.io
import bokeh.io.state
from bokeh.models import HoverTool, Panel, Tabs
from bokeh.plotting import figure, ColumnDataSource
from matplotlib import pyplot as plt
import numpy as np

from six.moves import cPickle
import scipy.stats

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import filepath
from six.moves import zip


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

			tfList = ["basal (no TF)"] + sorted(sim_data.tf_to_active_inactive_conditions)
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
			promoter_copy_number = rna_synth_prob_reader.readColumn("promoter_copy_number")

			tf_idx = tf_ids.index(tf)
			tf_targets = sim_data.tf_to_fold_change[tf]
			tf_target_indexes = np.array([
				rna_ids.index(tf_target + "[c]") for tf_target in tf_targets])

			tfTargetBoundCountsAll = n_bound_TF_per_TU[:, tf_target_indexes, tf_idx]
			tfTargetSynthProbAll = rna_synth_prob_reader.readColumn("rnaSynthProb")[:, tf_target_indexes]
			tf_target_promoter_copies_all = promoter_copy_number[:, tf_target_indexes]

			for i, tfTarget in enumerate(sorted(sim_data.tf_to_fold_change[tf])):
				tfTargetBoundCounts = tfTargetBoundCountsAll[:, i].reshape(-1)
				tf_target_copies = tf_target_promoter_copies_all[:, i].reshape(-1)

				expectedProbBound.append(sim_data.pPromoterBound[tf + "__" + tfStatus][tf])
				simulatedProbBound.append(
					(tfTargetBoundCounts[5:].astype(np.float64)/tf_target_copies[5:]).mean())

				tfTargetSynthProb = tfTargetSynthProbAll[:, i].reshape(-1)

				target_idx = rna_ids.index(tfTarget + "[c]")
				regulating_tf_idxs = delta_prob['deltaJ'][delta_prob['deltaI'] == target_idx]

				for j in regulating_tf_idxs:
					if tfTarget not in targetToTfType:
						targetToTfType[tfTarget] = []
					targetToTfType[tfTarget].append(sim_data.process.transcription_regulation.tf_to_tf_type[tf_ids[j]])

				expectedSynthProb.append(sim_data.process.transcription.rna_synth_prob[tf + "__" + tfStatus][target_idx])
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

		# Tools for all plots

		# Probability bound - hover for ID
		source1 = ColumnDataSource(data=dict(
			x=np.log10(expectedProbBound),
			y=np.log10(simulatedProbBound),
			ID=targetId,
			condition=targetCondition,
			))
		hover1 = HoverTool(tooltips=[("ID", "@ID"), ("condition", "@condition")])
		tools1 = [hover1, 'box_zoom', 'lasso_select', 'pan', 'wheel_zoom', 'undo', 'redo', 'reset']
		s1 = figure(
			x_axis_label="log10(Expected probability bound)",
			y_axis_label="log10(Simulated probability bound)",
			width=800,
			height=800,
			tools=tools1,
			)
		s1.scatter("x", "y", source=source1)

		html_dir = filepath.makedirs(plotOutDir, "html_plots")
		bokeh.io.output_file(os.path.join(html_dir, plotOutFileName + "__probBound" + ".html"), title=plotOutFileName)
		bokeh.io.save(s1)

		# Synthesis probability - hover for ID
		source2 = ColumnDataSource(data=dict(
			x=np.log10(expectedSynthProb),
			y=np.log10(simulatedSynthProb),
			ID=targetId,
			condition=targetCondition,
			))
		hover2 = HoverTool(tooltips=[("ID", "@ID"), ("condition", "@condition")])
		tools2 = [hover2, 'box_zoom', 'lasso_select', 'pan', 'wheel_zoom', 'undo', 'redo', 'reset']
		s2 = figure(
			x_axis_label="log10(Expected synthesis probability)",
			y_axis_label="log10(Simulated synthesis probability)",
			width=800,
			height=800,
			tools=tools2,
			)
		s2.scatter("x", "y", source=source2)

		bokeh.io.output_file(os.path.join(html_dir, plotOutFileName + "__synthProb" + ".html"), title=plotOutFileName)
		bokeh.io.save(s2)

		# Synthesis probability - filter targets by TF type
		bokeh.io.output_file(os.path.join(html_dir, plotOutFileName + "__synthProb__interactive" + ".html"), title=plotOutFileName)

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
		x0[np.where(tfTypes != "0CS")] = np.nan
		x1 = np.copy(expectedSynthProb)
		x1[np.where(tfTypes != "1CS")] = np.nan
		x2 = np.copy(expectedSynthProb)
		x2[np.where(tfTypes != "2CS")] = np.nan
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

		axis_min = np.floor(min(
			np.log10(expectedSynthProb)[np.isfinite(np.log10(expectedSynthProb))].min(),
			np.log10(simulatedSynthProb)[np.isfinite(np.log10(simulatedSynthProb))].min()))
		axis_max = np.ceil(max(
			np.log10(expectedSynthProb).max(),
			np.log10(simulatedSynthProb).max()))
		hover3 = HoverTool(tooltips=[("ID", "@ID"), ("condition", "@condition")])
		tools3 = [hover3, 'box_zoom', 'lasso_select', 'pan', 'wheel_zoom', 'undo', 'redo', 'reset']

		tabs = []
		data = [
			(expectedSynthProb, simulatedSynthProb, 'All'),
			(x0, y0, '0CS'),
			(x1, y1, '1CS'),
			(x2, y2, '2CS'),
			(x01, y01, '0CS and 1CS'),
			(x02, y02, '0CS and 2CS'),
			(x12, y12, '1CS and 2CS'),
			]
		for x, y, title in data:
			fig = figure(
				x_axis_label="log10(Expected synthesis probability)",
				y_axis_label="log10(Simulated synthesis probability)",
				plot_width=800,
				plot_height=800,
				x_range=(axis_min, axis_max),
				y_range=(axis_min, axis_max),
				tools=tools3,
				)
			source = ColumnDataSource(data=dict(
				x=np.log10(x),
				y=np.log10(y),
				ID=targetId,
				condition=targetCondition,
				))
			fig.scatter('x', 'y', source=source)
			tabs.append(Panel(child=fig, title=title))

		tab_plot = Tabs(tabs=tabs)
		bokeh.io.save(tab_plot)
		bokeh.io.state.curstate().reset()


if __name__ == "__main__":
	Plot().cli()
