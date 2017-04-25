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
	fig.set_figheight(6)

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

		if varIdx == 0:
			idxs_to_keep = np.where((0.65 < scaled_initial_masses) & (scaled_initial_masses < 1.25))[0]
		if varIdx == 1:
			idxs_to_keep = np.where((0.8 < scaled_initial_masses) & (scaled_initial_masses < 1.3))[0]
		if varIdx == 2:
			idxs_to_keep = np.where((0.6 < scaled_initial_masses) & (scaled_initial_masses < 1.2))[0]

		scaled_initial_masses = scaled_initial_masses[idxs_to_keep]
		scaled_added_masses = scaled_added_masses[idxs_to_keep]

		ax0 = plt.subplot2grid((1,3), (0,plotIdx))
		ax0.plot(scaled_initial_masses, scaled_added_masses, '.', color = "black", alpha = 0.2, zorder=1, markeredgewidth = 0.0)
		
		if varIdx == 0:
			nbins = 5
		elif varIdx == 1:
			nbins = 6
		if varIdx == 2:
			nbins = 20

		n_cell_cutoff = 5
		n, _ = np.histogram(scaled_initial_masses, bins=nbins)
		sy, _ = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses)
		sy2, _ = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses*scaled_added_masses)
		mean = sy / n
		std = np.sqrt(sy2/(n-1) - n*mean*mean/(n-1))
		ax0.errorbar(((_[1:] + _[:-1])/2)[n > n_cell_cutoff], mean[n > n_cell_cutoff], yerr=std[n > n_cell_cutoff], color = "black", linewidth=1, zorder=2)



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
		# ax0.text(0.6, 0.41, r"$m_{add}$=%.3f$\times$$m_{init}$ + %.3f"%(slope,intercept), fontsize=FONT_SIZE-2)
		# ax0.text(0.6, 0.35, "p-value={}".format(p_value), fontsize=FONT_SIZE-2)
		factor = 1.58*std_err
		#ax0.text(0.6, 0.4, "Slope 95 percent CI: %.3f - %.3f"%(slope - factor, slope+factor))

		ax0.set_title(
			title_list[varIdx] + ", n=%d, n*=%d" % ((len(all_cells) - fail), len(scaled_initial_masses)) + "\n" +
			r"$m_{add}$=%.3f$\times$$m_{init}$ + %.3f" % (slope,intercept) + "\n" +
			"p-value=%0.2g" % p_value,
			fontsize=FONT_SIZE)

		print "Slope 95 percent CI: %.3f -> %.3f"%(slope - factor, slope+factor)
		print "Mean: {}".format(n[n>n_cell_cutoff].mean())
		# ax0.errorbar(sj_mean_x, sj_mean_y, sj_error)

		#ax0.axhline(1., linewidth = 1, color = "black", alpha = 0.9)
		ax0.set_ylim([0.45, 1.5])
		# ax0.set_xlim([0.6, 1.4])

		ax0.get_yaxis().get_major_formatter().set_useOffset(False)
		ax0.get_xaxis().get_major_formatter().set_useOffset(False)

		# ax0.set_title("n = {}".format(n_cells))
		if varIdx == 1:
			ax0.set_ylabel("Normed added mass", fontsize=FONT_SIZE)
		ax0.set_xlabel("Normed initial mass", fontsize=FONT_SIZE)

		plt.subplots_adjust(bottom = 0.2)

		whitePadSparklineAxis(ax0)

		for tick in ax0.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		for tick in ax0.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)


	####### STRIPPED FIGURES #########

	title_list = ["Glucose minimal\n" + r"$\tau = $" + "44 min", "Glucose minimal anaerobic\n" + r"$\tau = $" + "100 min", "Glucose minimal + 20 amino acids\n" + r"$\tau = $" + "22 min"]

	for varIdx in [0,1,2]:
		fig = plt.figure()
		fig.set_figwidth(1.73)
		fig.set_figheight(1.18)

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

		if varIdx == 0:
			idxs_to_keep = np.where((0.65 < scaled_initial_masses) & (scaled_initial_masses < 1.25))[0]
		if varIdx == 1:
			idxs_to_keep = np.where((0.8 < scaled_initial_masses) & (scaled_initial_masses < 1.3))[0]
		if varIdx == 2:
			idxs_to_keep = np.where((0.6 < scaled_initial_masses) & (scaled_initial_masses < 1.2))[0]

		scaled_initial_masses = scaled_initial_masses[idxs_to_keep]
		scaled_added_masses = scaled_added_masses[idxs_to_keep]

		ax0 = plt.subplot2grid((1,1), (0,0))
		ax0.plot(scaled_initial_masses, scaled_added_masses, '.', color = "black", alpha = 0.2, zorder=1, markeredgewidth = 0.0)
		
		if varIdx == 0:
			nbins = 5
		elif varIdx == 1:
			nbins = 6
		if varIdx == 2:
			nbins = 20

		n_cell_cutoff = 5
		n, _ = np.histogram(scaled_initial_masses, bins=nbins)
		sy, _ = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses)
		sy2, _ = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses*scaled_added_masses)
		mean = sy / n
		std = np.sqrt(sy2/(n-1) - n*mean*mean/(n-1))
		ax0.errorbar(((_[1:] + _[:-1])/2)[n > n_cell_cutoff], mean[n > n_cell_cutoff], yerr=std[n > n_cell_cutoff], color = "black", linewidth=1, zorder=2)

		ax0.set_title(title_list[varIdx] + ", n=%d, n*=%d"% (len(all_cells) - fail, len(scaled_initial_masses)), fontsize=FONT_SIZE)

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
		#ax0.text(0.6, 0.41, r"$m_{add}$=%.3f$\times$$m_{init}$ + %.3f"%(slope,intercept), fontsize=FONT_SIZE-2)
		#ax0.text(0.6, 0.35, "p-value={}".format(p_value), fontsize=FONT_SIZE-2)
		factor = 1.58*std_err
		#ax0.text(0.6, 0.4, "Slope 95 percent CI: %.3f - %.3f"%(slope - factor, slope+factor))

		print "Slope 95 percent CI: %.3f -> %.3f"%(slope - factor, slope+factor)
		print "Mean: {}".format(n[n>n_cell_cutoff].mean())
		# ax0.errorbar(sj_mean_x, sj_mean_y, sj_error)

		#ax0.axhline(1., linewidth = 1, color = "black", alpha = 0.9)
		ax0.set_ylim([0.45, 1.5])
		# ax0.set_xlim([0.6, 1.4])

		ax0.get_yaxis().get_major_formatter().set_useOffset(False)
		ax0.get_xaxis().get_major_formatter().set_useOffset(False)

		# ax0.set_title("n = {}".format(n_cells))
		if varIdx == 1:
			ax0.set_ylabel("Normed added mass", fontsize=FONT_SIZE)
		ax0.set_xlabel("Normed initial mass", fontsize=FONT_SIZE)

		plt.subplots_adjust(bottom = 0.2)

		whitePadSparklineAxis(ax0)

		for tick in ax0.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 
		for tick in ax0.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE) 

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
		exportFigure(plt, plotOutDir, plotOutFileName + str(varIdx) + "_stripped" ,metadata, transparent = True)
		plt.close("all")

	# x-axis histograms
	fig = plt.figure()
	fig.set_figwidth(11.)
	fig.set_figheight(6)

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

		ax0 = plt.subplot2grid((1,3), (0,varIdx))
		ax0.hist(scaled_initial_masses, 50)

		if varIdx == 0:
			ax0.axvline(x = 0.65, color = "r", linestyle = "--")
			ax0.axvline(x = 1.25, color = "r", linestyle = "--")
			ax0.set_title(title_list[varIdx] + "\n" + "[0.65, 1.25]", fontsize = FONT_SIZE)


		if varIdx == 1:
			ax0.axvline(x = 0.8, color = "r", linestyle = "--")
			ax0.axvline(x = 1.3, color = "r", linestyle = "--")
			ax0.set_title(title_list[varIdx] + "\n" + "[0.8, 1.3]", fontsize = FONT_SIZE)

		if varIdx == 2:
			ax0.axvline(x = 0.6, color = "r", linestyle = "--")
			ax0.axvline(x = 1.2, color = "r", linestyle = "--")
			ax0.set_title(title_list[varIdx] + "\n" + "[0.6, 1.2]", fontsize = FONT_SIZE)

		ax0.set_xlabel("Normed initial mass", fontsize=FONT_SIZE)

		plt.subplots_adjust(bottom = 0.2)

		whitePadSparklineAxis(ax0)

		for tick in ax0.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in ax0.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_histogram_scaled_initial_mass" ,metadata, transparent = True)
	plt.close("all")

	# y-axis histograms
	fig = plt.figure()
	fig.set_figwidth(11)
	fig.set_figheight(6)

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
		ax0.hist(scaled_added_masses, 50)

		if varIdx == 0:
			ax0.axvline(x = 0.45, color = "r", linestyle = "--")
			ax0.axvline(x = 1.5, color = "r", linestyle = "--")
			ax0.set_title(title_list[varIdx] + "\n" + "[0.45, 1.5]", fontsize = FONT_SIZE)

		if varIdx == 1:
			ax0.axvline(x = 0.45, color = "r", linestyle = "--")
			ax0.axvline(x = 1.5, color = "r", linestyle = "--")
			ax0.set_title(title_list[varIdx] + "\n" + "[0.45, 1.5]", fontsize = FONT_SIZE)

		if varIdx == 2:
			ax0.axvline(x = 0.45, color = "r", linestyle = "--")
			ax0.axvline(x = 1.5, color = "r", linestyle = "--")
			ax0.set_title(title_list[varIdx] + "\n" + "[0.45, 1.5]", fontsize = FONT_SIZE)

		ax0.set_xlabel("Normed added mass", fontsize=FONT_SIZE)

		plt.subplots_adjust(bottom = 0.2)

		whitePadSparklineAxis(ax0)

		for tick in ax0.yaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)
		for tick in ax0.xaxis.get_major_ticks():
			tick.label.set_fontsize(FONT_SIZE)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName + "_histogram_scaled_added_mass" ,metadata, transparent = True)
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
