#!/usr/bin/env python
"""
Plots fluxes across different conditions.

@author: Sam Bray
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/19/2017'
"""
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

	if not os.path.isdir(inputDir):
		raise Exception, "inputDir does not currently exist as a directory"

	ap = AnalysisPaths(inputDir, variant_plot = True)
	variants = sorted(ap._path_data['variant'].tolist()) # Sorry for accessing private data

	
	if len(variants) <= 1:
		return

	all_cells = sorted(ap.get_cells(variant = variants, seed = [0], generation = [0]))

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	#make structures to hold mean flux values
	mean_fluxes = []
	BURN_IN_STEPS = 20
	n_variants=0
	IDs=[]

	#Puts you into the specific simulation's data.  Pull fluxes from here  #TODO LEARN HOW TO PULL FLUXES FROM LISTENER FILE (see kineticsflux comparison)
	for variant, simDir in zip(variants, all_cells):
		sim_data = cPickle.load(open(ap.get_variant_kb(variant), "rb"))
		simOutDir = os.path.join(simDir, "simOut")

		#crafting area
		enzymeKineticsReader = TableReader(os.path.join(simOutDir, "FBAResults"))# "EnzymeKinetics"))
		actualFluxes = enzymeKineticsReader.readColumn("reactionFluxes")#"actualFluxes")
		IDs = enzymeKineticsReader.readAttribute("reactionIDs")
		enzymeKineticsReader.close()

		actualAve = np.mean(actualFluxes[BURN_IN_STEPS:, :], axis = 0)
		mean_fluxes.append(actualAve)
		n_variants=n_variants+1
		

	###Plot the fluxes
	plt.figure(figsize = (8.5, 11))

	#Generalizred plotting
	for j in range (0,n_variants):
		for k in range (0,n_variants):
			if j<=k:
				continue
			plt.subplot(n_variants-1, n_variants-1, j+k)
			plt.plot(np.log10(mean_fluxes[j][:]),np.log10(mean_fluxes[k][:]),'o')
			plt.plot([-12, 0], [-12, 0], color='k', linestyle='-', linewidth=2)
			plt.xlabel('Variant ' + str(j) + ' Flux')
			plt.ylabel('Variant ' + str(k) + ' Flux')
			plt.ylim((-11,0))
			plt.xlim((-11,0))


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


	#nifty fun tool
	# Bokeh
	if len(mean_fluxes) < 2:
		return

	from bokeh.plotting import figure, output_file, ColumnDataSource, show
	from bokeh.models import HoverTool, BoxZoomTool, LassoSelectTool, PanTool, WheelZoomTool, ResizeTool, UndoTool, RedoTool

	# Plot first metabolite to initialize plot settings
	x = np.log10(mean_fluxes[0][:])
	y = np.log10(mean_fluxes[1][:])
	

	source = ColumnDataSource(
		data = dict(
			x = x,
			y = y,
			rxn = IDs)
		)

	hover = HoverTool(
		tooltips = [
			("ID", "@rxn"),
			]
		)
	
	TOOLS = [hover, 
		BoxZoomTool(),
		LassoSelectTool(),
		PanTool(),
		WheelZoomTool(),
		ResizeTool(),
		UndoTool(),
		RedoTool(),
		 "reset"
		 ]

	p = figure(x_axis_label = "Variant 0 Flux", 
		y_axis_label = "Variant 1 Flux",
		width = 800,
		height = 800,
		tools = TOOLS,
		)

	p.circle('x', 'y', size=5, source=source)#np.log10(mean_fluxes[0][:]),np.log10(mean_fluxes[1][:]), size=10)
	p.line([-12,0], [-12,0], color="firebrick", line_width=2)


	if not os.path.exists(os.path.join(plotOutDir, "html_plots")):
		os.makedirs(os.path.join(plotOutDir, "html_plots"))

	import bokeh.io
	bokeh.io.output_file(os.path.join(plotOutDir, "html_plots", plotOutFileName + ".html"), title=plotOutFileName, autosave=False)
	bokeh.io.save(p)
	bokeh.io.curstate().reset()


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
