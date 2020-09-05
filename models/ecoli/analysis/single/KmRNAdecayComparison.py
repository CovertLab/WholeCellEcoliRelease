"""
Plots counts of rna degraded and the resulting free NMPs

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/2015 - Updated 8/10/2015
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle

from wholecell.analysis.analysis_tools import exportFigure
from wholecell.analysis.analysis_tools import read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))

		if sim_data.constants.EndoRNaseCooperation:
			KmFirstOrderDecay = sim_data.process.rna_decay.Km_first_order_decay
			KmNonLinearDecay = (sim_data.process.transcription.rna_data['Km_endoRNase'].asNumber())

			# Compute deviation
			Error = np.average(np.abs(KmFirstOrderDecay
								- KmNonLinearDecay)
								/ KmFirstOrderDecay
					* 100)

			# Plotting
			plt.figure(figsize = (6, 12))

			plt.subplot(3,1,1)
			plt.loglog(KmFirstOrderDecay, KmNonLinearDecay, 'o', markeredgecolor = 'k', markerfacecolor = 'none')

			minLine = np.round(1.2 * min((np.log10(KmFirstOrderDecay)).min(), (np.log10(KmNonLinearDecay).min())))
			plt.loglog([np.power(10, minLine), 1], [np.power(10, minLine), 1], '--r')

			plt.xlabel("Km First Order Decay (Log10, M)", fontsize = 14)
			plt.ylabel("Km Non-linear Decay (Log10, M)", fontsize = 14)
			plt.title("Relative error = %.2f%%" % Error, fontsize = 16)
			# print(np.corrcoef(KmFirstOrderDecay, KmNonLinearDecay)[0,1])


		if sim_data.constants.EndoRNaseCooperation:
			plt.subplot(3,1,2)
			GprimeKm = sim_data.process.rna_decay.Km_convergence
			FprimeKm = np.log10(1 - GprimeKm[GprimeKm < 1])
			plt.hist(FprimeKm)

			plt.ylabel("Number of genes", fontsize = 14)
			plt.xlabel("Log10(1 - g\'(Km))", fontsize = 14)
			PercentageConvergence = len(GprimeKm[GprimeKm < 1.]) / float(len(GprimeKm)) * 100.
			plt.title("Convergence of %.0f%% Km\'s" % PercentageConvergence, fontsize = 16)


		# Sensitivity analysis kcatEndoRNases
		# TODO: does this ever get set and should it be a variant analysis plot?
		if sim_data.constants.SensitivityAnalysisKcatEndo:
			cellDensity = sim_data.constants.cellDensity
			cellVolume = sim_data.mass.avg_cell_dry_mass_init / cellDensity / sim_data.mass.cell_dry_mass_fraction
			countsToMolar = 1 / (sim_data.constants.n_Avogadro * cellVolume)

			rnaIds = sim_data.process.transcription.rna_data["id"]
			(rna_counts_bulk,) = read_bulk_molecule_counts(simOutDir, rnaIds)
			RNAcounts = rna_counts_bulk[-1, :]

			ax = plt.subplot(3,1,3)
			ax.set_xscale("log", nonposx='clip')
			ax.set_yscale("log", nonposy='clip')

			ConvergenceFactor = 10

			fractionRNAkm_avg = []
			fractionRNAkm_sd = []
			Kcats = []
			Km = []

			for kcat in sim_data.process.rna_decay.sensitivity_analysis_kcat:
				KMcounts = 1 / countsToMolar.asNumber() * sim_data.process.rna_decay.sensitivity_analysis_kcat[kcat]
				ResIni = sim_data.process.rna_decay.sensitivity_analysis_kcat_res_ini[kcat]
				ResOpt = sim_data.process.rna_decay.sensitivity_analysis_kcat_res_opt[kcat]

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

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
