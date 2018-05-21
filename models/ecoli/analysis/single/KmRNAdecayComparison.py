#!/usr/bin/env python
"""
Plots counts of rna degraded and the resulting free NMPs

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/2015 - Updated 8/10/2015
"""

from __future__ import division

import argparse
import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

import wholecell.utils.constants
from wholecell.io.tablereader import TableReader

import scipy.stats as st

FONT = {
		'size'	:	14
		}


def setAxisMaxMin(axis, data):
	ymax = np.max(data)
	ymin = 0
	if ymin == ymax:
		axis.set_yticks([ymin])
	else:
		axis.set_yticks([ymin, ymax])

def sparklineAxis(axis, x, y, tickPos, lineType, color):
	axis.plot(x, y, linestyle = 'steps' + lineType, color = color, linewidth = 2)
	axis.spines['top'].set_visible(False)
	axis.spines['bottom'].set_visible(False)
	axis.yaxis.set_ticks_position(tickPos)
	axis.xaxis.set_ticks_position('none')
	axis.tick_params(which = 'both', direction = 'out')
	axis.tick_params(labelbottom = 'off')
	for tl in axis.get_yticklabels():
		tl.set_color(color)


def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))

	if sim_data.constants.EndoRNaseCooperation:
		KmFirstOrderDecay = sim_data.process.rna_decay.KmFirstOrderDecay
		KmNonLinearDecay = (sim_data.process.transcription.rnaData["KmEndoRNase"].asNumber())

		FC = np.log10(1 - (KmNonLinearDecay / KmFirstOrderDecay))

		# Compute deviation
		Error = np.average(np.abs(KmFirstOrderDecay
							- KmNonLinearDecay)
							/ KmFirstOrderDecay
				* 100)

		# Plotting
		plt.figure(figsize = (6, 12))
		matplotlib.rc('font', **FONT)
		max_yticks = 5


		plt.subplot(3,1,1)
		plt.loglog(KmFirstOrderDecay, KmNonLinearDecay, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

		minLine = np.round(1.2 * min((np.log10(KmFirstOrderDecay)).min(), (np.log10(KmNonLinearDecay).min())))
		plt.loglog([np.power(10, minLine), 1], [np.power(10, minLine), 1], '--r')

		plt.xlabel("Km First Order Decay (Log10, M)", fontsize = 14)
		plt.ylabel("Km Non-linear Decay (Log10, M)", fontsize = 14)
		plt.title("Relative error = %.2f%%" % Error, fontsize = 16)
		print np.corrcoef(KmFirstOrderDecay, KmNonLinearDecay)[0,1]


	if sim_data.constants.EndoRNaseCooperation:
		plt.subplot(3,1,2)
		GprimeKm = sim_data.process.rna_decay.KmConvergence
		FprimeKm = np.log10(1 - GprimeKm[GprimeKm < 1])
		plt.hist(FprimeKm)

		plt.ylabel("Number of genes", fontsize = 14)
		plt.xlabel("Log10(1 - g\'(Km))", fontsize = 14)
		PercentageConvergence = len(GprimeKm[GprimeKm < 1.]) / float(len(GprimeKm)) * 100.
		plt.title("Convergence of %.0f%% Km\'s" % PercentageConvergence, fontsize = 16)



	# Sensitivity analysis kcatEndoRNases
	cellDensity = sim_data.constants.cellDensity
	cellVolume = sim_data.mass.avgCellDryMassInit / cellDensity / sim_data.mass.cellDryMassFraction
	countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

	isMRna = sim_data.process.transcription.rnaData["isMRna"]
	rnaIds = sim_data.process.transcription.rnaData["id"]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)
	rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]
	bulkMolecules.close()
	RNAcounts = rnaCountsBulk[-1, :]

	if sim_data.constants.SensitivityAnalysisKcatEndo:
		ax = plt.subplot(3,1,3)
		ax.set_xscale("log", nonposx='clip')
		ax.set_yscale("log", nonposy='clip')

		width = 1
		ConvergenceFactor = 10


		fractionRNAkm_avg = []
		fractionRNAkm_sd = []
		Kcats = []
		Km = []

		for kcat in sim_data.process.rna_decay.SensitivityAnalysisKcat:
			KMcounts = 1 / countsToMolar.asNumber() * sim_data.process.rna_decay.SensitivityAnalysisKcat[kcat]
			ResIni = sim_data.process.rna_decay.SensitivityAnalysisKcat_ResIni[kcat]
			ResOpt = sim_data.process.rna_decay.SensitivityAnalysisKcat_ResOpt[kcat]

			if ResIni > ResOpt * ConvergenceFactor:
				Kcats = np.append(Kcats, kcat)
				Km = np.append(Km, np.average(KMcounts))
				fractionRNAkm_avg = np.append(fractionRNAkm_avg, np.average(RNAcounts / KMcounts))
				fractionRNAkm_sd = np.append(fractionRNAkm_sd, np.std(RNAcounts / KMcounts) )

		plt.errorbar(Kcats, fractionRNAkm_avg, yerr = [np.zeros(len(fractionRNAkm_sd)), fractionRNAkm_sd], fmt = 'o')

		z = np.polyfit(np.log10(Kcats), np.log10(fractionRNAkm_avg), 1)
		p = np.poly1d(z)

		minX = np.log10(Kcats.min() / 2)
		maxX = np.log10(np.round(2 * Kcats.max()))
		minY = np.power(10, p(minX))
		maxY = np.power(10, p(maxX))
		plt.plot([np.power(10, minX), np.power(10,maxX)], [minY, maxY], '--k')

		plt.xlabel("Kcat EndoRNase (1/s)", fontsize = 14)
		plt.ylabel("RNA / Km", fontsize = 14)
		plt.title("Endo-nucleolytic cleavage operating on linear regime \n(Michaelis-Menten law, RNA/Km << 1)", fontsize = 12)


	plt.subplots_adjust(hspace = 0.4, wspace = 0.2)

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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "None")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
