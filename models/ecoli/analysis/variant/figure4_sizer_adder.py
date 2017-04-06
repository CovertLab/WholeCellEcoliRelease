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

from scipy.stats import linregress
def mm2inch(value):
	return value * 0.0393701

FONT_SIZE=9
trim = 0.05

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return r_value**2

def main(inputDir, plotOutDir, plotOutFileName, validationDataFile = None, metadata = None):
	return
	if not os.path.isdir(inputDir):
		raise Exception, "variantDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(inputDir, variant_plot = True)

	if ap.n_generation == 1:
		print "Need more data to create addedMass"
		return

	fig = plt.figure()
	fig.set_figwidth(11)
	fig.set_figheight(4)

	title_list = ["Glucose minimal\n" + r"$\tau = $" + "44 min", "Glucose minimal anaerobic\n" + r"$\tau = $" + "100 min", "Glucose minimal + 20 amino acids\n" + r"$\tau = $" + "22 min"]

	for varIdx in [0,1,2]:

		if varIdx == 0:
			plotIdx = 1
		elif varIdx == 1:
			plotIdx = 0
		elif varIdx == 2:
			plotIdx = 2

		initial_masses = np.zeros(0)
		final_masses = np.zeros(0)

		if varIdx == 0:
			gen = [2,3]
		elif varIdx == 1:
			gen = [2,3]
		elif varIdx == 2:
			gen = [6,7]

		all_cells = ap.get_cells(generation=gen, variant=[varIdx])
		fail = 0
		for simDir in all_cells:
			try:
				simOutDir = os.path.join(simDir, "simOut")
				mass = TableReader(os.path.join(simOutDir, "Mass"))
				cellMass = mass.readColumn("dryMass")

				initial_masses = np.hstack((initial_masses, cellMass[0]))
				final_masses = np.hstack((final_masses, cellMass[-1]))
			except:
				fail+=1

		print final_masses.size
		added_masses = final_masses - initial_masses

		scaled_initial_masses = initial_masses / initial_masses.mean()
		scaled_added_masses = added_masses / added_masses.mean()

		ax0 = plt.subplot2grid((1,3), (0,plotIdx))
		ax0.plot(scaled_initial_masses, scaled_added_masses, '.', color = "grey", alpha = 0.5, zorder=1)
		
		nbins = 5
		if varIdx == 2:
			nbins = 10

		n_cell_cutoff = 5
		n, _ = np.histogram(scaled_initial_masses, bins=nbins)
		sy, _ = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses)
		sy2, _ = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses*scaled_added_masses)
		mean = sy / n
		std = np.sqrt(sy2/(n-1) - n*mean*mean/(n-1))
		ax0.errorbar(((_[1:] + _[:-1])/2)[n > n_cell_cutoff], mean[n > n_cell_cutoff], yerr=std[n > n_cell_cutoff], color = "black", linewidth=1, zorder=2)

		ax0.set_title(title_list[varIdx] + ", n={}".format(len(all_cells) - fail), fontsize=FONT_SIZE)

		# z = np.polyfit(scaled_initial_masses, scaled_added_masses, 1)
		# p = np.poly1d(z)
		# ax0.plot(scaled_initial_masses, p(scaled_initial_masses), color = "blue")
		# ax0.text(0.6, 1.5, r"$m_{add}$=%.3f$\times$$m_{init}$ + %.3f"%(z[0],z[1]))

		# x = scaled_initial_masses
		# y = scaled_added_masses
		# yhat = p(x)                         # or [p(z) for z in x]
		# ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
		# ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
		# sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
		# r_squared = ssreg / sstot

		slope, intercept, r_value, p_value, std_err = linregress(scaled_initial_masses, scaled_added_masses)
		ax0.plot(scaled_initial_masses, slope * scaled_initial_masses + intercept, color = "blue")
		ax0.text(0.6, 0.41, r"$m_{add}$=%.3f$\times$$m_{init}$ + %.3f"%(slope,intercept), fontsize=FONT_SIZE-2)
		ax0.text(0.6, 0.35, "p-value={}".format(p_value), fontsize=FONT_SIZE-2)
		factor = 1.58*std_err
		#ax0.text(0.6, 0.4, "Slope 95 percent CI: %.3f - %.3f"%(slope - factor, slope+factor))

		print "Slope 95 percent CI: %.3f -> %.3f"%(slope - factor, slope+factor)
		print "Mean: {}".format(n[n>n_cell_cutoff].mean())
		# ax0.errorbar(sj_mean_x, sj_mean_y, sj_error)

		#ax0.axhline(1., linewidth = 1, color = "black", alpha = 0.9)
		ax0.set_ylim([0.3, 1.7])
		ax0.set_xlim([0.6, 1.4])

		ax0.get_yaxis().get_major_formatter().set_useOffset(False)
		ax0.get_xaxis().get_major_formatter().set_useOffset(False)

		# ax0.set_title("n = {}".format(n_cells))
		ax0.set_ylabel("Normed added mass", fontsize=FONT_SIZE)
		ax0.set_xlabel("Normed initial mass", fontsize=FONT_SIZE)

		plt.subplots_adjust(bottom = 0.2, wspace= 0.6)

		whitePadSparklineAxis(ax0)

		for tick in ax0.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		for tick in ax0.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)



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
