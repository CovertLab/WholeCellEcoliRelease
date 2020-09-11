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

		if sim_data.constants.endoRNase_cooperation:
			width = 1

			LossKm = sim_data.process.rna_decay.stats_fit['LossKm']
			LossKmOpt = sim_data.process.rna_decay.stats_fit['LossKmOpt']
			RnegKmOpt = sim_data.process.rna_decay.stats_fit['RnegKmOpt']
			ResKm = sim_data.process.rna_decay.stats_fit['ResKm']
			ResKmOpt = sim_data.process.rna_decay.stats_fit['ResKmOpt']
			ResEndoRNKm = sim_data.process.rna_decay.stats_fit['ResEndoRNKm']
			ResEndoRNKmOpt = sim_data.process.rna_decay.stats_fit['ResEndoRNKmOpt']
			ResScaledKm = sim_data.process.rna_decay.stats_fit['ResScaledKm']
			ResScaledKmOpt = sim_data.process.rna_decay.stats_fit['ResScaledKmOpt']

			StatsFit = []
			ScoreNames = []

			# Sensitivity analysis alpha
			Residuals = sim_data.process.rna_decay.sensitivity_analysis_alpha_residual

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

			rnaDegRates = sim_data.process.transcription.rna_data['deg_rate']

			cellDensity = sim_data.constants.cell_density
			cellVolume = sim_data.mass.avg_cell_dry_mass_init / cellDensity / sim_data.mass.cell_dry_mass_fraction
			countsToMolar = 1 / (sim_data.constants.n_avogadro * cellVolume)

			rnaIds = sim_data.process.transcription.rna_data["id"]
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
