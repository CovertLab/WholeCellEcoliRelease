from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot
from wholecell.io.tablereader import TableReader

from wholecell.utils.sparkline import whitePadSparklineAxis

from scipy.stats import linregress


FONT_SIZE = 5

INIT_MASS_LOWER_LIM = 0.6
INIT_MASS_UPPER_LIM = 1.6
ADDED_MASS_LOWER_LIM = 0
ADDED_MASS_UPPER_LIM = 1.7

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot = True)

		if ap.n_generation == 1:
			print "Need more data to create addedMass"
			return

		allScatter = plt.figure()
		allScatter.set_figwidth(11)
		allScatter.set_figheight(6)

		plt.style.use('seaborn-deep')
		color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

		title_list = [r"Glucose minimal, $\tau = $44 min", r"Glucose minimal anaerobic, $\tau = $100 min", r"Glucose minimal + 20 amino acids, $\tau = $25 min"]

		for varIdx in ap.get_variants():

			if varIdx == 0:
				plotIdx = 1
				gen = [2,3]
			elif varIdx == 1:
				plotIdx = 0
				gen = [2,3]
			elif varIdx == 2:
				plotIdx = 2
				gen = [2,3]
			else:
				continue

			initial_masses = np.zeros(0)
			final_masses = np.zeros(0)

			all_cells = ap.get_cells(generation=gen, variant=[varIdx])
			if len(all_cells) == 0:
				continue

			fail = 0
			for simDir in all_cells:
				try:
					simOutDir = os.path.join(simDir, "simOut")
					mass = TableReader(os.path.join(simOutDir, "Mass"))
					cellMass = mass.readColumn("dryMass")

					initial_masses = np.hstack((initial_masses, cellMass[0]))
					final_masses = np.hstack((final_masses, cellMass[-1]))
				except Exception as e:
					print e
					fail+=1

			added_masses = final_masses - initial_masses

			scaled_initial_masses = initial_masses / initial_masses.mean()
			scaled_added_masses = added_masses / added_masses.mean()

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
				title_list[varIdx] + ", n=%d" % ((len(all_cells) - fail), ) + "\n" +
				r"$m_{add}$=%.3f$\times$$m_{init}$ + %.3f" % (slope,intercept) + "\n" +
				"r-value=%0.2g" % r_value + "\n" +
				"p-value=%0.2g" % p_value,
				fontsize=FONT_SIZE)

			ax.set_xlim([INIT_MASS_LOWER_LIM, INIT_MASS_UPPER_LIM])
			ax.set_ylim([ADDED_MASS_LOWER_LIM, ADDED_MASS_UPPER_LIM])
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
			fig.set_figwidth(3)
			fig.set_figheight(2)
			ax = plt.subplot2grid((1,1), (0,0))
			ax.plot(scaled_initial_masses, scaled_added_masses, '.', color = color_cycle[0], alpha = 0.25, ms=6, zorder=1, markeredgewidth = 0.0, clip_on=False)
			ax.plot(scaled_initial_masses, slope * scaled_initial_masses + intercept, color = 'k')

			ax.set_xlim([INIT_MASS_LOWER_LIM, INIT_MASS_UPPER_LIM])
			ax.set_ylim([ADDED_MASS_LOWER_LIM, ADDED_MASS_UPPER_LIM])

			ax.get_yaxis().get_major_formatter().set_useOffset(False)
			ax.get_xaxis().get_major_formatter().set_useOffset(False)

			whitePadSparklineAxis(ax)

			ax.tick_params(which='both', bottom=True, left=True,
				top=False, right=False, labelbottom=True, labelleft=True,
				labelsize=FONT_SIZE)

			ax.set_xlabel("")
			ax.set_ylabel("")

			plt.tight_layout()
			exportFigure(plt, plotOutDir, plotOutFileName + str(varIdx) + "_stripped", metadata, transparent = True)

		plt.figure(allScatter.number)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
