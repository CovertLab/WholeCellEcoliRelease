"""
Plot to assess sensitivity of pabB behavior to model parameters.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/1/17
"""

from __future__ import absolute_import, division, print_function

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.subgen_expression import FACTORS
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.utils.sparkline import whitePadSparklineAxis


CONC_UNITS = units.mol / units.L
THRESHOLD = 0.999  # 99.9% of target conc
ENZYME_IDS = ['PABASYN-CPLX[c]', 'PABSYNMULTI-CPLX[c]']  # both catalyze reaction of interest
METABOLITE_ID = 'METHYLENE-THF[c]'

FONTSIZE = 10
LABELSIZE = 8
MARKERSIZE = 1


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, "inputDir does not currently exist as a directory"
		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get cells
		ap = AnalysisPaths(inputDir, variant_plot = True)
		if ap.n_variant != len(FACTORS):
			print("This plot expects all variants of subgen_expression")
			return

		# Get constants from wildtype variant
		sim_data = cPickle.load(open(ap.get_variant_kb(4), "rb")) # 4 is the wildtype variant
		cellDensity = sim_data.constants.cellDensity
		nAvogadro = sim_data.constants.nAvogadro
		metabolite_target = sim_data.process.metabolism.concDict[METABOLITE_ID]
		metabolite_threshold = (THRESHOLD * metabolite_target).asNumber(CONC_UNITS)

		# Investigate each variant
		enzyme_depletion = np.zeros([ap.n_seed, ap.n_variant])
		metabolite_depletion = np.zeros([ap.n_seed, ap.n_variant])

		for variant in xrange(ap.n_variant):
			for seed in xrange(ap.n_seed):
				cells = ap.get_cells(variant=[variant], seed=[seed])
				time_enzyme_depleted = []  # seconds
				time_metabolite_depleted = []  # seconds

				for i, simDir in enumerate(cells):
					simOutDir = os.path.join(simDir, "simOut")

					main_reader = TableReader(os.path.join(simOutDir, "Main"))
					mass_reader = TableReader(os.path.join(simOutDir, "Mass"))

					# Get molecule counts
					(enzyme_counts, metabolite_counts) = read_bulk_molecule_counts(simOutDir, (ENZYME_IDS, [METABOLITE_ID]))

					# Compute time with zero counts of enzyme
					time_step_sec = main_reader.readColumn("timeStepSec")
					time_enzyme_depleted.append(time_step_sec[np.sum(enzyme_counts, axis=1) == 0].sum())

					# Compute time with end products under the target concentration
					mass = units.fg * mass_reader.readColumn("cellMass")
					volume = mass / cellDensity
					metabolite_conc = (1 / nAvogadro / volume * metabolite_counts).asNumber(CONC_UNITS)
					time_metabolite_depleted.append(time_step_sec[metabolite_conc < metabolite_threshold].sum())

				# Record MENE-CPLX depletion
				total_time = main_reader.readColumn("time")[-1] + time_step_sec[-1]
				fraction_enzyme_depleted = np.sum(time_enzyme_depleted) / total_time
				enzyme_depletion[seed, variant] = fraction_enzyme_depleted

				# Record end product depletion
				fraction_metabolite_depleted = np.sum(time_metabolite_depleted) / total_time
				metabolite_depletion[seed, variant] = fraction_metabolite_depleted

		# Compute average and standard deviations
		metabolite_depletion_avg = np.average(metabolite_depletion, axis = 0)
		metabolite_depletion_std = np.std(metabolite_depletion, axis = 0)
		enzyme_depletion_avg = np.average(enzyme_depletion, axis = 0)
		enzyme_depletion_std = np.std(enzyme_depletion, axis = 0)

		# Plot
		fig, axesList = plt.subplots(2, 1, figsize = (8, 8))
		ax1, ax2 = axesList
		xvals = np.arange(ap.n_variant)
		fig.suptitle("Sensitivity Analysis: pabB depletion")

		for ax, avg, std in zip(axesList, [metabolite_depletion_avg, enzyme_depletion_avg], [metabolite_depletion_std, enzyme_depletion_std]):
			ax.scatter(xvals, avg, edgecolor = "none", clip_on = False, s = MARKERSIZE)
			ax.errorbar(xvals, avg, yerr = std, color = "b", linewidth = 1, clip_on = False, fmt = "o", capsize = 4, capthick = 1, markeredgecolor = "none")

		ax1.set_title("Enzyme depletion", fontsize = FONTSIZE)
		ax2.set_title("Metabolite depletion", fontsize = FONTSIZE)
		xlabels = ["1/10 x", "1/8 x", "1/4 x", "1/2 x", "1 x", "2 x", "4 x", "8 x", "10 x"]
		title_tags = ["counts = 0", "<%s%% of wildtype" % (THRESHOLD * 100)]
		for i, ax in enumerate([ax1, ax2]):
			ax.set_ylabel("Fraction of Time\n%s" % title_tags[i], fontsize = FONTSIZE)
			ax.set_xlabel("Factor of change of pabB synthesis probability", fontsize = FONTSIZE)
			ax.set_xlim([-0.25, 8.25])
			whitePadSparklineAxis(ax)
			ax.set_xticks(xvals)
			ax.set_xticklabels(xlabels)
			ax.set_yticks([0, 1])

		plt.subplots_adjust(hspace = 1, wspace = 1, top = 0.9, bottom = 0.1)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")

		# Plot clean versions for figure
		FIRST = True
		for avg, std, filename in zip([metabolite_depletion_avg, enzyme_depletion_avg], [metabolite_depletion_std, enzyme_depletion_std], ["pabB", "methylene-thf"]):
			fig, ax = plt.subplots(1, 1, figsize = (10, 3))
			ax.scatter(xvals, avg, edgecolor = "none", clip_on = False, s = MARKERSIZE)
			ax.errorbar(xvals, avg, yerr = std, color = "b", linewidth = 1, clip_on = False, fmt = "o", capsize = 4, capthick = 1, markeredgecolor = "none")
			ax.set_xlim([-0.25, 8.25])
			if FIRST:
				FIRST = False
				whitePadSparklineAxis(ax, False)
			else:
				whitePadSparklineAxis(ax)
			ax.set_xticks(xvals)
			ax.set_xticklabels([])
			ax.set_yticks([0, 1])
			ax.set_yticklabels([])
			exportFigure(plt, plotOutDir, plotOutFileName + "_%s" % filename, metadata)
			plt.close("all")


if __name__ == "__main__":
	Plot().cli()
