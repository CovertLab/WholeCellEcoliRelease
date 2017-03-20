#!/usr/bin/env python

import argparse
import os
import re
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from wholecell.utils.sparkline import whitePadSparklineAxis

def mm2inch(value):
	return value * 0.0393701

FONT_SIZE=9
trim = 0.05

def main(variantDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile = None, metadata = None):
	if not os.path.isdir(variantDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get all cells in each seed
	ap = AnalysisPaths(variantDir, cohort_plot = True)

	if ap.n_generation == 1:
		print "Need more data to create addedMass"
		return

	initial_masses = np.zeros(0)
	final_masses = np.zeros(0)

	all_cells = ap.get_cells(generation=[4,5,6,7])

	if all_cells.tolist() == []:
		return

	for simDir in all_cells:
		simOutDir = os.path.join(simDir, "simOut")
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = mass.readColumn("dryMass")

		initial_masses = np.hstack((initial_masses, cellMass[0]))
		final_masses = np.hstack((final_masses, cellMass[-1]))

	print final_masses.size
	added_masses = final_masses - initial_masses

	scaled_initial_masses = initial_masses / initial_masses.mean()
	scaled_added_masses = added_masses / added_masses.mean()


	# # Sucjoon data
	# sj_mean = [(0.599573863268619, 1.0424971716015852),
	# 			(0.7025146209428339, 1.013345898239791),
	# 			(0.7957469172398974, 1.0228576014663655),
	# 			(0.898666963529841, 1.0065888091212933),
	# 			(1.0015818319737173, 0.9935406370304014),
	# 			(1.0980399263710288, 0.9966059219025472),
	# 			(1.2009340834306335, 0.9964402308283771),
	# 			(1.2909664708577877, 0.9962952511384784),
	# 			(1.397076070325505, 0.9961243822182404),
	# 			(1.499913271078364, 1.0313855139400554),
	# 			(1.5931611009136308, 1.0312353564040888)]

	# sj_mean_x = [x[0] for x in sj_mean]
	# sj_mean_y = [y[1] for y in sj_mean]

	# sj_error_top = [(0.599320148811296, 1.2003075640564278),
	# 				(0.7023489298686638, 1.1164057463735655),
	# 				(0.7988277356502467, 1.1065885502289898),
	# 				(0.9017684933244621, 1.0774372768671956),
	# 				(1.0046730060762024, 1.0708303452846646),
	# 				(1.0979260137575368, 1.0674595674945173),
	# 				(1.1975995505629613, 1.0705196745205958),
	# 				(1.290821491167889, 1.0864726182555315),
	# 				(1.3969000235591993, 1.1056254708603763),
	# 				(1.499747580004194, 1.1344453620738304),
	# 				(1.5929695206091217, 1.1503983058087661)]

	# sj_error_bottom = [(0.5998327555720095, 0.8814661588925619),
	# 					(0.7059009322711844, 0.9070602520057676),
	# 					(0.7991125171839766, 0.9294544362490642),
	# 					(0.8987912318354687, 0.929293923020962),
	# 					(1.0016853888950736, 0.929128231946792),
	# 					(1.0981641946766563, 0.9193110358022161),
	# 					(1.2010738852744645, 0.9094834839655048),
	# 					(1.2943372486479348, 0.8996714656669966),
	# 					(1.397246939245743, 0.8898439138302852),
	# 					(1.5000737843064662, 0.9315462860604613),
	# 					(1.5933578590642077, 0.9088517867452315)]

	# sj_error =[sj_error_top[idx][1] - sj_error_bottom[idx][1] for idx in range(len(sj_mean))]

	#### FIGURE 4D ####
	fig, ax0 = plt.subplots(1)
	mult = 2
	fig.set_figwidth(mm2inch(40)*mult)
	fig.set_figheight(mm2inch(38)*mult)

	ax0.plot(scaled_initial_masses, scaled_added_masses, '.', color = "grey", alpha = 0.5, zorder=1)
	nbins = 30
	n_cell_cutoff = 10
	n, _ = np.histogram(scaled_initial_masses, bins=nbins)
	sy, _ = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses)
	sy2, _ = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses*scaled_added_masses)
	mean = sy / n
	std = np.sqrt(sy2/(n-1) - n*mean*mean/(n-1))
	ax0.errorbar(((_[1:] + _[:-1])/2)[n > n_cell_cutoff], mean[n > n_cell_cutoff], yerr=std[n > n_cell_cutoff], color = "black", linewidth=2, zorder=2)

	# ax0.errorbar(sj_mean_x, sj_mean_y, sj_error)

	ax0.axhline(1., linewidth = 1, color = "black", alpha = 0.9)
	ax0.text(np.max(ax0.get_xlim()) - 0.001, 1., "adder")
	ax0.set_ylim([0., 2.])
	ax0.set_xlim([0.6, 1.4])

	ax0.get_yaxis().get_major_formatter().set_useOffset(False)
	ax0.get_xaxis().get_major_formatter().set_useOffset(False)

	# ax0.set_title("n = {}".format(n_cells))
	ax0.set_ylabel("Scaled added mass " + r"$(\frac{m_{added}}{\langle m_{added} \rangle})$", fontsize=FONT_SIZE)
	ax0.set_xlabel("Scaled initial mass " + r"$(\frac{m_{initial}}{\langle m_{initial} \rangle})$", fontsize=FONT_SIZE)

	plt.subplots_adjust(left = 0.2, bottom = 0.2, wspace= 0.6)

	whitePadSparklineAxis(ax0)

	for tick in ax0.yaxis.get_major_ticks():
		tick.label.set_fontsize(FONT_SIZE) 
	for tick in ax0.xaxis.get_major_ticks():
		tick.label.set_fontsize(FONT_SIZE) 

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_f", metadata)

	for axes in [ax0]:
		axes.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom='off',      # ticks along the bottom edge are off
			top='off',         # ticks along the top edge are off
			labelbottom='off') # labels along the bottom edge are off
		axes.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			left='off',      # ticks along the bottom edge are off
			right='off',         # ticks along the top edge are off
			labelleft='off') # labels along the bottom edge are off

		axes.set_xlabel("")
		axes.set_ylabel("")

	plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_f_stripped" ,metadata, transparent = True)
	plt.close("all")

	#### FIGURE 4E ####

	fig, ax1 = plt.subplots(1)
	mult = 2
	fig.set_figwidth(mm2inch(40)*mult)
	fig.set_figheight(mm2inch(38)*mult)

	# Plot contours for all but first generation
	# H, xedges, yedges = np.histogram2d(initial_masses, added_masses, bins=np.round(n_cells/10))

	H, xedges, yedges = np.histogram2d(initial_masses, added_masses, bins=5)

	X, Y = np.meshgrid(xedges, yedges)
	ax1.contour(X[:-1,:-1], Y[:-1,:-1], H.transpose(), cmap="Greys")
	ax1.get_yaxis().get_major_formatter().set_useOffset(False)

	ax1.set_xlabel("Initial mass (pg)", fontsize=FONT_SIZE)
	ax1.set_ylabel("Added mass (pg)", fontsize=FONT_SIZE)
	whitePadSparklineAxis(ax1)
	plt.subplots_adjust(left = 0.2, bottom = 0.2, wspace= 0.6)


	for tick in ax1.yaxis.get_major_ticks():
		tick.label.set_fontsize(FONT_SIZE) 
	for tick in ax1.xaxis.get_major_ticks():
		tick.label.set_fontsize(FONT_SIZE) 

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_g", metadata)

	for axes in [ax1]:
		axes.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom='off',      # ticks along the bottom edge are off
			top='off',         # ticks along the top edge are off
			labelbottom='off') # labels along the bottom edge are off
		axes.tick_params(
			axis='y',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			left='off',      # ticks along the bottom edge are off
			right='off',         # ticks along the top edge are off
			labelleft='off') # labels along the bottom edge are off

		axes.set_xlabel("")
		axes.set_ylabel("")

	plt.subplots_adjust(top = 1, bottom = trim, left = trim, right = 1)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_g_stripped" ,metadata, transparent = True)
	plt.close("all")

	#### UNUSED PLOT FOR NOW ####
	# Linear mapping

	# ax2.plot(initial_masses, final_masses, '.', color = "black")
	# z = np.polyfit(initial_masses, final_masses, 1)
	# p = np.poly1d(z)
	# ax2.plot(initial_masses, p(initial_masses), '--', color = "grey")
	# text_x = np.min(ax2.get_xlim())
	# text_y = np.max(ax2.get_ylim())
	# ax2.text(text_x, text_y, r"$m_f$=%.3f$\times$$m_i$ + %.3f"%(z[0],z[1]))

	# ax2.set_xlabel("Initial mass (pg)")
	# ax2.set_ylabel("Final mass (pg)")

	# whitePadSparklineAxis(ax2)
	# plt.subplots_adjust(left = 0.2, bottom = 0.2, wspace= 0.6)

	# from wholecell.analysis.analysis_tools import exportFigure
	# exportFigure(plt, plotOutDir, plotOutFileName + "_e", metadata)
	# plt.close("all")



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
