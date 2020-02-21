"""
Plots counts of rna degraded and the resulting free NMPs

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/2015 - Updated 8/10/2015
"""

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from models.ecoli.analysis import singleAnalysisPlot

FONT = {
		'size'	:	10
		}


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))

		ntp_ids = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
		endoRnaseIds = sim_data.process.rna_decay.endoRnaseIds
		exoRnaseIds = sim_data.moleculeGroups.exoRnaseIds
		RnaseIds = np.concatenate((endoRnaseIds, exoRnaseIds))

		(RnaseCounts, exoRnaseCounts, endoRnaseCounts, ntpCounts) = read_bulk_molecule_counts(
			simOutDir, (RnaseIds, exoRnaseIds, endoRnaseIds, ntp_ids))

		rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
		FractionActiveEndoRNases = rnaDegradationListenerFile.readColumn('FractionActiveEndoRNases')
		DiffRelativeFirstOrderDecay = rnaDegradationListenerFile.readColumn('DiffRelativeFirstOrderDecay')
		fragmentBasesDigested = rnaDegradationListenerFile.readColumn('fragmentBasesDigested')
		rnaDegradationListenerFile.close()

		TranscriptElongationListenerFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
		countNTPsUSed = TranscriptElongationListenerFile.readColumn('countNTPsUSed')
		countRnaSynthesized = TranscriptElongationListenerFile.readColumn('countRnaSynthesized')
		TranscriptElongationListenerFile.close()

		totalexoRnaseCounts = exoRnaseCounts.sum(axis = 1)
		totalendoRnaseCounts = endoRnaseCounts.sum(axis = 1)

		# Load metabolism production
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime
		deltaMetabolites = fbaResults.readColumn("deltaMetabolites")
		outputMoleculeIDs = np.array(fbaResults.readAttribute("metaboliteNames"))
		fbaResults.close()

		# Plotting
		plt.figure(figsize = (8.5, 11))
		plt.rc('font', **FONT)
		max_yticks = 5
		n_rows = 6
		n_cols = 2

		ax = plt.subplot(n_rows, n_cols, 1)
		plt.plot(time / 60., countRnaSynthesized.sum(axis = 1))
		plt.ylabel("RNAs synthesized", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(n_rows, n_cols, 3)
		plt.plot(time / 60., countRnaDegraded.sum(axis = 1))
		plt.ylabel("RNAs degraded", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(n_rows, n_cols, 5)
		plt.plot(time / 60., totalendoRnaseCounts)
		plt.ylabel("EndoRNase counts", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(n_rows, n_cols, 2)
		plt.plot(time / 60., countNTPsUSed / 1e6)
		plt.ylabel("Transcription ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("NTPs needed(x$10^{%d}$) = %.2f" % (6, (countNTPsUSed.sum() / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(n_rows, n_cols, 7)
		plt.plot(time / 60., totalexoRnaseCounts)
		plt.ylabel("ExoRNase counts", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		IdxAtp = (np.where("ATP[c]" == outputMoleculeIDs))[0][0]; ATP = np.sum(deltaMetabolites[:, IdxAtp])
		IdxGtp = (np.where("GTP[c]" == outputMoleculeIDs))[0][0]; GTP = np.sum(deltaMetabolites[:, IdxGtp])
		IdxCtp = (np.where("CTP[c]" == outputMoleculeIDs))[0][0]; CTP = np.sum(deltaMetabolites[:, IdxCtp])
		IdxUtp = (np.where("UTP[c]" == outputMoleculeIDs))[0][0]; UTP = np.sum(deltaMetabolites[:, IdxUtp])
		NtpsProduced = ATP + GTP + CTP + UTP
		ax = plt.subplot(n_rows, n_cols, 4)
		plt.plot(time / 60., (deltaMetabolites[:, IdxAtp] + deltaMetabolites[:, IdxGtp] + deltaMetabolites[:, IdxCtp] + deltaMetabolites[:, IdxUtp]) / 1e6)
		plt.ylabel("Metabolism ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("NTPs produced (x$10^{%d}$) = %.2f" % (6, (sum(deltaMetabolites[:, IdxAtp] + deltaMetabolites[:, IdxGtp] + deltaMetabolites[:, IdxCtp] + deltaMetabolites[:, IdxUtp])  / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(n_rows, n_cols, 9)
		plt.plot(time / 60., FractionActiveEndoRNases * 100)
		plt.ylabel("EndoRN capacity (%)", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(n_rows, n_cols, 6)
		plt.plot(time / 60., fragmentBasesDigested / 1e6)
		plt.ylabel("Exo-digestion ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("NTPs recycled (x$10^{%d}$) = %.2f" % (6, (fragmentBasesDigested.sum() / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(n_rows, n_cols, 12)
		plt.plot(time / 60., DiffRelativeFirstOrderDecay)
		plt.xlabel("Time (min)")
		plt.ylabel("sum(Residuals)", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(n_rows, n_cols, 8)
		plt.plot(time / 60., (ntpCounts[:, 0] + ntpCounts[:, 1] + ntpCounts[:, 2] + ntpCounts[:, 3]) / 1e6)
		plt.ylabel("Net production ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("NTPs required for cell division (x$10^{%d}$) = %.2f" % (6, ((ntpCounts[0, 0] + ntpCounts[0, 1] + ntpCounts[0, 2] + ntpCounts[0, 3])  / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		# compute active ExoRNase capacity (%)
		ActiveExoRNcapacity = fragmentBasesDigested.astype(float) / (totalexoRnaseCounts * sim_data.constants.KcatExoRNase.asNumber()) * 100

		ax = plt.subplot(n_rows, n_cols, 11)
		plt.plot(time / 60., ActiveExoRNcapacity)
		plt.xlabel("Time (min)")
		plt.ylabel("ExoRN capacity (%)", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		# compute instantaneous balance of nTPs
		InstantaneousNTPs = - countNTPsUSed + (deltaMetabolites[:, IdxAtp] + deltaMetabolites[:, IdxGtp] + deltaMetabolites[:, IdxCtp] + deltaMetabolites[:, IdxUtp]) + fragmentBasesDigested

		ax = plt.subplot(n_rows, n_cols, 10)
		plt.plot(time / 60., InstantaneousNTPs / 1e6)
		plt.ylabel("Balance ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("Average instantaneous balance (x$10^{%d}$) = %.4f" % (6, (np.mean(InstantaneousNTPs) / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)


		plt.subplots_adjust(hspace = 0.6, wspace = 0.35)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
