from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.io.tablereader import TableReader

from wholecell.utils.sparkline import whitePadSparklineAxis

from scipy.stats import linregress


FONT_SIZE=9
trim = 0.05


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot = True)

		if ap.n_generation == 1:
			print("Need more data to create addedMass")
			return

		allScatter = plt.figure()
		allScatter.set_figwidth(11)
		allScatter.set_figheight(6)

		xHist = plt.figure()
		xHist.set_figwidth(11)
		xHist.set_figheight(6)

		yHist = plt.figure()
		yHist.set_figwidth(11)
		yHist.set_figheight(6)

		title_list = ["Glucose minimal\n" + r"$\tau = $" + "44 min", "Glucose minimal anaerobic\n" + r"$\tau = $" + "100 min", "Glucose minimal + 20 amino acids\n" + r"$\tau = $" + "22 min"]

		plot = False
		for varIdx in ap.get_variants():

			if varIdx == 0:
				plotIdx = 1
				gen = [2,3]
			elif varIdx == 1:
				plotIdx = 0
				gen = [2,3]
			elif varIdx == 2:
				plotIdx = 2
				gen = [6,7]
			else:
				continue

			initial_masses = np.zeros(0)
			final_masses = np.zeros(0)

			all_cells = ap.get_cells(generation=gen, variant=[varIdx])
			if len(all_cells) == 0:
				continue
			plot = True

			fail = 0
			for simDir in all_cells:
				try:
					simOutDir = os.path.join(simDir, "simOut")
					mass = TableReader(os.path.join(simOutDir, "Mass"))
					cellMass = mass.readColumn("dryMass")

					initial_masses = np.hstack((initial_masses, cellMass[0]))
					final_masses = np.hstack((final_masses, cellMass[-1]))
				except Exception as e:
					print(e)
					fail+=1

			added_masses = final_masses - initial_masses

			all_scaled_initial_masses = initial_masses / initial_masses.mean()
			all_scaled_added_masses = added_masses / added_masses.mean()

			idxs_to_keep = np.where((0.6 < all_scaled_initial_masses) & (all_scaled_initial_masses < 1.25) & (0.45 < all_scaled_added_masses) & (all_scaled_added_masses < 1.5))

			scaled_initial_masses = all_scaled_initial_masses[idxs_to_keep]
			scaled_added_masses = all_scaled_added_masses[idxs_to_keep]

			nbins = 5

			n, xbin = np.histogram(scaled_initial_masses, bins=nbins)
			sy, xbin = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses)
			sy2, xbin = np.histogram(scaled_initial_masses, bins=nbins, weights=scaled_added_masses*scaled_added_masses)
			mean = sy / n
			std = np.sqrt(sy2/(n-1) - n*mean*mean/(n-1))

			slope, intercept, r_value, p_value, std_err = linregress(scaled_initial_masses, scaled_added_masses)

			# plot all scatter plots
			plt.figure(allScatter.number)
			ax = plt.subplot2grid((1,3), (0,plotIdx))
			ax.plot(scaled_initial_masses, scaled_added_masses, '.', color = "black", alpha = 0.2, zorder=1, markeredgewidth = 0.0)
			ax.errorbar(((xbin[1:] + xbin[:-1])/2), mean, yerr=std, color = "black", linewidth=1, zorder=2)
			ax.plot(scaled_initial_masses, slope * scaled_initial_masses + intercept, color = "blue")

			ax.set_title(
				title_list[varIdx] + ", n=%d, n*=%d" % ((len(all_cells) - fail), len(scaled_initial_masses)) + "\n" +
				r"$m_{add}$=%.3f$\times$$m_{init}$ + %.3f" % (slope,intercept) + "\n" +
				"p-value=%0.2g" % p_value,
				fontsize=FONT_SIZE)

			ax.set_xlim([0.6, 1.25])
			ax.set_ylim([0.45, 1.5])
			ax.get_yaxis().get_major_formatter().set_useOffset(False)
			ax.get_xaxis().get_major_formatter().set_useOffset(False)

			if varIdx == 1:
				ax.set_ylabel("Normed added mass", fontsize=FONT_SIZE)
			ax.set_xlabel("Normed initial mass", fontsize=FONT_SIZE)

			plt.subplots_adjust(bottom = 0.2)

			whitePadSparklineAxis(ax)

			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)
			for tick in ax.xaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)

			# plot stripped figure
			fig = plt.figure()
			fig.set_figwidth(1.73)
			fig.set_figheight(1.18)
			ax = plt.subplot2grid((1,1), (0,0))
			ax.plot(scaled_initial_masses, scaled_added_masses, '.', color = "black", alpha = 0.2, zorder=1, markeredgewidth = 0.0)
			ax.errorbar(((xbin[1:] + xbin[:-1])/2), mean, yerr=std, color = "black", linewidth=1, zorder=2)
			ax.set_title(title_list[varIdx] + ", n=%d, n*=%d"% (len(all_cells) - fail, len(scaled_initial_masses)), fontsize=FONT_SIZE)
			ax.plot(scaled_initial_masses, slope * scaled_initial_masses + intercept, color = "blue")

			ax.set_ylim([0.45, 1.5])

			ax.get_yaxis().get_major_formatter().set_useOffset(False)
			ax.get_xaxis().get_major_formatter().set_useOffset(False)

			plt.subplots_adjust(bottom = 0.2)

			whitePadSparklineAxis(ax)

			ax.tick_params(
				axis='x',
				which='both',
				bottom=False,
				top=False,
				labelbottom=False)
			ax.tick_params(
				axis='y',
				which='both',
				left=False,
				right=False,
				labelleft=False)

			ax.set_xlabel("")
			ax.set_ylabel("")

			plt.subplots_adjust(top = 0.95, bottom = 3 * trim, left = 2 * trim, right = 0.95, hspace = 0, wspace = 0)

			exportFigure(plt, plotOutDir, plotOutFileName + str(varIdx) + "_stripped", metadata, transparent = True)

			# plot histogram for x-axis
			plt.figure(xHist.number)
			bins = 50
			ax = plt.subplot2grid((1,3), (0,plotIdx))
			ax.hist(all_scaled_initial_masses, bins)

			ax.axvline(x = 0.6, color = "r", linestyle = "--")
			ax.axvline(x = 1.25, color = "r", linestyle = "--")
			ax.set_title(title_list[varIdx] + "\n" + "[0.6, 1.25]", fontsize = FONT_SIZE)

			ax.set_xlabel("Normed initial mass", fontsize = FONT_SIZE)

			plt.subplots_adjust(bottom = 0.2)

			whitePadSparklineAxis(ax)

			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)
			for tick in ax.xaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)

			# plot histogram for y-axis
			plt.figure(yHist.number)
			ax = plt.subplot2grid((1,3), (0,plotIdx))
			ax.hist(all_scaled_added_masses, bins)

			ax.axvline(x = 0.45, color = "r", linestyle = "--")
			ax.axvline(x = 1.5, color = "r", linestyle = "--")
			ax.set_title(title_list[varIdx] + "\n" + "[0.45, 1.5]", fontsize = FONT_SIZE)

			ax.set_xlabel("Normed added mass", fontsize = FONT_SIZE)

			plt.subplots_adjust(bottom = 0.2)

			whitePadSparklineAxis(ax)

			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)
			for tick in ax.xaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)

		if plot:
			plt.figure(allScatter.number)
			exportFigure(plt, plotOutDir, plotOutFileName, metadata)
			plt.figure(xHist.number)
			exportFigure(plt, plotOutDir, plotOutFileName + "_histogram_scaled_initial_mass" ,metadata, transparent = True)
			plt.figure(yHist.number)
			exportFigure(plt, plotOutDir, plotOutFileName + "_histogram_scaled_added_mass" ,metadata, transparent = True)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
