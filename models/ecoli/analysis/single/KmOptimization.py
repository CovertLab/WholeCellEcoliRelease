#!/usr/bin/env python
"""
@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/27/2016
"""

from __future__ import division

import argparse
import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

import cPickle

THRESHOLD = 1e-13 # roughly, the mass of an electron

FG_PER_DALTON = 1.6605402e-9

# TODO: get these from the KB
REPRESENTATIVE_MASSES = {
	"proton":1.007 * FG_PER_DALTON,
	"amino acid":109 * FG_PER_DALTON,
	"ATP":551 * FG_PER_DALTON,
	"protein":40e3 * FG_PER_DALTON,
	"ribosome":2700e3 * FG_PER_DALTON
	}

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	if sim_data.constants.EndoRNaseCooperation:

		mass = TableReader(os.path.join(simOutDir, "Mass"))

		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		processMassDifferences = mass.readColumn("processMassDifferences")

		processNames = mass.readAttribute("processNames")

		mass.close()

		avgProcessMassDifferences = np.abs(processMassDifferences).sum(axis = 0) / len(time)

		index = np.arange(len(processNames))

		width = 1

		sim_data = cPickle.load(open(simDataFile, "rb"))

		LossKm = sim_data.process.rna_decay.StatsFit['LossKm']
		LossKmOpt = sim_data.process.rna_decay.StatsFit['LossKmOpt']
		RnegKmOpt = sim_data.process.rna_decay.StatsFit['RnegKmOpt']
		ResKm = sim_data.process.rna_decay.StatsFit['ResKm']
		ResKmOpt = sim_data.process.rna_decay.StatsFit['ResKmOpt']
		ResEndoRNKm = sim_data.process.rna_decay.StatsFit['ResEndoRNKm']
		ResEndoRNKmOpt = sim_data.process.rna_decay.StatsFit['ResEndoRNKmOpt']
		ResScaledKm = sim_data.process.rna_decay.StatsFit['ResScaledKm']
		ResScaledKmOpt = sim_data.process.rna_decay.StatsFit['ResScaledKmOpt']

		StatsFit = []
		ScoreNames = []

		# Sensitivity analysis alpha
		Residuals = sim_data.process.rna_decay.SensitivityAnalysisAlphaResidual

		for alpha in sorted(Residuals):
			ScoreNames.append('Residuals rescaled, alpha = ' + str(alpha))
			StatsFit.append(Residuals[alpha])

		StatsFit = StatsFit + [0, ResScaledKmOpt, ResScaledKm, ResEndoRNKmOpt, ResEndoRNKm, ResKmOpt, ResKm, RnegKmOpt, LossKmOpt, LossKm]
		ScoreNames = ScoreNames + [
			'',
			'Residuals rescaled(KmOpt), M/s',
			'Residuals rescaled(Km), M/s',
			'Residuals EndoRN(KmOpt)',
			'Residuals EndoRN(Km)',
			'Residuals(KmOpt)',
			'Residuals(Km)',
			'Total Negative Regularization(KmOpt)',
			'Total Loss(KmOpt)',
			'Total Loss(Km)',
			]

		index = np.arange(len(StatsFit))

		plt.figure(figsize = (8.5, 11))

		axes = plt.axes()

		r1 = axes.barh(index, StatsFit, width, log = True, color = (0.2, 0.2, 0.9))

		axes.set_yticks(index+width/2)
		axes.set_yticklabels(ScoreNames) #, rotation = -45)

		# If a THRESHOLD is defined:
		# axes.plot([THRESHOLD, THRESHOLD], [index[0], index[-1]+width], 'k--', linewidth=3)
		# plt.text(THRESHOLD, index[-1], "electron", rotation = "vertical", va = "center", ha = "right")

		rnaDegRates = sim_data.process.transcription.rnaData['degRate']

		cellDensity = sim_data.constants.cellDensity
		cellVolume = sim_data.mass.avgCellDryMassInit / cellDensity / sim_data.mass.cellDryMassFraction
		countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

		rnaIds = sim_data.process.transcription.rnaData["id"]
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)
		rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]
		rnaCountsInitial = rnaCountsBulk[-1, :]
		rnaConcInitial = countsToMolar * rnaCountsInitial
		rnaDecay = rnaConcInitial * rnaDegRates

		REPRESENTATIVE_SCORES = {
				'Sum (Kd * RNAcounts), M/s': np.sum(rnaDecay).asNumber(),
				#'Sum (Kd), 1/s': np.sum(rnaDegRates.asNumber()),
			}

		for name in REPRESENTATIVE_SCORES:
			mass = REPRESENTATIVE_SCORES[name]
			plt.axvline(mass, color = "k", linestyle = "dashed", linewidth = "3")
			plt.text(mass, index[-1], name, rotation = "vertical", ha = "right")

		plt.title("Scores Km's non-linear optimization")

		plt.tight_layout()
		plt.grid(True, which = "major")

		from wholecell.analysis.analysis_tools import exportFigure
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
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
	parser.add_argument("--validationDataFile", help = "KB file name", type = str, default = "None")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
