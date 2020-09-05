"""
Plots frequency of observing at least 1 transcript during a cell's life.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/10/2017
"""

from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle

import bokeh.io
import bokeh.io.state
from bokeh.models import HoverTool
from bokeh.plotting import figure, ColumnDataSource
import matplotlib.pyplot as plt
import numpy as np

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.analysis import multigenAnalysisPlot
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import filepath


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		# Get mRNA data
		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnaIds = sim_data.process.transcription.rna_data["id"]
		isMRna = sim_data.process.transcription.rna_data['is_mRNA']
		mRnaIds = np.where(isMRna)[0]
		mRnaNames = np.array([rnaIds[x] for x in mRnaIds])

		# Get whether or not mRNAs were transcribed
		transcribedBool = []
		simulatedSynthProbs = []
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")

			rnaSynthProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
			simulatedSynthProb = np.mean(rnaSynthProb.readColumn("rnaSynthProb")[:, mRnaIds], axis = 0)
			rnaSynthProb.close()
			simulatedSynthProbs.append(simulatedSynthProb)

			mRNA_counts_reader = TableReader(
				os.path.join(simOutDir, 'mRNACounts'))
			moleculeCounts = mRNA_counts_reader.readColumn("mRNA_counts")
			moleculeCountsSumOverTime = moleculeCounts.sum(axis = 0)
			mRnasTranscribed = np.array([x != 0 for x in moleculeCountsSumOverTime])
			transcribedBool.append(mRnasTranscribed)

		transcribedBool = np.array(transcribedBool)
		simulatedSynthProbs = np.array(simulatedSynthProbs)

		# Plot frequency vs. synthesis prob
		fig = plt.figure(figsize = (12, 12))
		ax = plt.subplot(1, 1, 1)
		# ax.scatter(np.log10(mRnaSynthProb), np.mean(transcribedBool, axis = 0))
		ax.scatter(np.log10(np.mean(simulatedSynthProbs, axis = 0)), np.mean(transcribedBool, axis = 0))
		ax.set_title("Frequency of observing at least 1 transcript in 1 generation", fontsize = 14)
		ax.set_xlabel("log10 (Transcript synthesis probability)", fontsize = 10)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		# Plot freq as histogram
		nGens = len(allDir)
		fig = plt.figure(figsize = (12, 12))
		ax = plt.subplot(1, 1, 1)
		ax.hist(np.mean(transcribedBool, axis = 0), nGens)
		ax.set_xlabel("Frequency of observing at least 1 transcript in 1 generation", fontsize = 14)
		ax.set_ylabel("Number of genes", fontsize = 14)
		exportFigure(plt, plotOutDir, plotOutFileName + "__histogram", metadata)
		plt.close("all")

		hover = HoverTool(tooltips = [("ID", "@ID")])
		tools = [hover, "box_zoom", "pan", "wheel_zoom", "tap", "save", "reset"]
		plot = figure(
			x_axis_label="log10 (transcript synthesis probability)",
			y_axis_label="Frequency of observing at least 1 transcript",
			width=800,
			height=800,
			tools=tools,
			)
		source = ColumnDataSource(data=dict(
			x=np.log10(np.mean(simulatedSynthProbs, axis=0)),
			y=np.mean(transcribedBool, axis=0),
			ID=mRnaNames,
			))
		plot.scatter("x", "y", source = source)

		html_dir = filepath.makedirs(plotOutDir, "html_plots")
		bokeh.io.output_file(os.path.join(html_dir, plotOutFileName + ".html"), title=plotOutFileName)
		bokeh.io.save(plot)
		bokeh.io.state.curstate().reset()


if __name__ == "__main__":
	Plot().cli()
