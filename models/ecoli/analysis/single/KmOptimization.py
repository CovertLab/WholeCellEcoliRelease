"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/27/2016
"""

from __future__ import absolute_import, division, print_function

from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))

		if sim_data.constants.EndoRNaseCooperation:
			width = 1

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
			axes.set_yticklabels(ScoreNames)

			rnaDegRates = sim_data.process.transcription.rnaData['degRate']

			cellDensity = sim_data.constants.cellDensity
			cellVolume = sim_data.mass.avgCellDryMassInit / cellDensity / sim_data.mass.cellDryMassFraction
			countsToMolar = 1 / (sim_data.constants.nAvogadro * cellVolume)

			rnaIds = sim_data.process.transcription.rnaData["id"]
			(rnaCountsBulk,) = read_bulk_molecule_counts(simOutDir, rnaIds)
			rnaCountsInitial = rnaCountsBulk[-1, :]
			rnaConcInitial = countsToMolar * rnaCountsInitial
			rnaDecay = rnaConcInitial * rnaDegRates

			REPRESENTATIVE_SCORES = {
					'Sum (Kd * RNAcounts), M/s': np.sum(rnaDecay).asNumber(),
					#'Sum (Kd), 1/s': np.sum(rnaDegRates.asNumber()),
				}

			for name in REPRESENTATIVE_SCORES:
				mass = REPRESENTATIVE_SCORES[name]
				plt.axvline(mass, color = "k", linestyle = "dashed", linewidth = 3)
				plt.text(mass, index[-1], name, rotation = "vertical", ha = "right")

			plt.title("Scores Km's non-linear optimization")

			plt.tight_layout()
			plt.grid(True, which = "major")

			exportFigure(plt, plotOutDir, plotOutFileName, metadata)
			plt.close("all")


if __name__ == "__main__":
	Plot().cli()
