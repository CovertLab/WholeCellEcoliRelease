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
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FONT = {
		'size'	:	10
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


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))

		endoRnaseIds = sim_data.process.rna_decay.endoRnaseIds
		exoRnaseIds = sim_data.moleculeGroups.exoRnaseIds
		RnaseIds = np.concatenate((endoRnaseIds, exoRnaseIds))

		# Load count data for s30 proteins, rRNA, and final 30S complex
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")

		# Get indexes
		proteinIndexes = np.array([moleculeIds.index(protein) for protein in RnaseIds], np.int)
		exoproteinIndexes = np.array([moleculeIds.index(protein) for protein in exoRnaseIds], np.int)
		endoproteinIndexes = np.array([moleculeIds.index(protein) for protein in endoRnaseIds], np.int)

		# Load data
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		RnaseCounts = bulkMolecules.readColumn("counts")[:, proteinIndexes]

		exoRnaseCounts = bulkMolecules.readColumn("counts")[:, exoproteinIndexes]
		endoRnaseCounts = bulkMolecules.readColumn("counts")[:, endoproteinIndexes]
		bulkMolecules.close()

		rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
		nucleotidesFromDegradation = rnaDegradationListenerFile.readColumn('nucleotidesFromDegradation')
		FractionActiveEndoRNases = rnaDegradationListenerFile.readColumn('FractionActiveEndoRNases')
		DiffRelativeFirstOrderDecay = rnaDegradationListenerFile.readColumn('DiffRelativeFirstOrderDecay')
		FractEndoRRnaCounts = rnaDegradationListenerFile.readColumn('FractEndoRRnaCounts')
		fragmentBasesDigested = rnaDegradationListenerFile.readColumn('fragmentBasesDigested')
		rnaDegradationListenerFile.close()

		TranscriptElongationListenerFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
		countNTPsUSed = TranscriptElongationListenerFile.readColumn('countNTPsUSed')
		countRnaSynthesized = TranscriptElongationListenerFile.readColumn('countRnaSynthesized')
		TranscriptElongationListenerFile.close()

		totalRnaseCounts = RnaseCounts.sum(axis = 1)
		requiredRnaseTurnover = nucleotidesFromDegradation / RnaseCounts.sum(axis = 1)

		totalexoRnaseCounts = exoRnaseCounts.sum(axis = 1)
		totalendoRnaseCounts = endoRnaseCounts.sum(axis = 1)

		# Load data
		growthLimitsDataFile = TableReader(os.path.join(simOutDir, "GrowthLimits"))

		# Translation
		gtpUsed = growthLimitsDataFile.readColumn("gtpAllocated")
		growthLimitsDataFile.close()

		# Load metabolism production
		fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		deltaMetabolites = fbaResults.readColumn("deltaMetabolites")
		outputMoleculeIDs = np.array(fbaResults.readAttribute("metaboliteNames"))
		fbaResults.close()

		# Load ntps required for cell doubling
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		NTP_IDS = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
		ntpIndexes = np.array([moleculeIds.index(ntpId) for ntpId in NTP_IDS], np.int)
		ntpCounts = bulkMolecules.readColumn("counts")[:, ntpIndexes]
		bulkMolecules.close()


		# Plotting
		plt.figure(figsize = (8.5, 11))
		plt.rc('font', **FONT)
		max_yticks = 5

		ax = plt.subplot(7,2,1)
		plt.plot(time / 60., countRnaSynthesized.sum(axis = 1))
		plt.ylabel("RNAs synthesized", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,2)
		plt.plot(time / 60., gtpUsed / 1e6)
		plt.ylabel("Translation ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("GTPs needed (x$10^{%d}$) = %.2f" % (6, (gtpUsed.sum() / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,3)
		plt.plot(time / 60., countRnaDegraded.sum(axis = 1))
		plt.ylabel("RNAs degraded", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,5)
		plt.plot(time / 60., totalendoRnaseCounts)
		plt.ylabel("EndoRNase counts", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,4)
		plt.plot(time / 60., countNTPsUSed / 1e6)
		plt.ylabel("Transcription ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("NTPs needed(x$10^{%d}$) = %.2f" % (6, (countNTPsUSed.sum() / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,7)
		plt.plot(time / 60., totalexoRnaseCounts)
		plt.ylabel("ExoRNase counts", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		IdxAtp = (np.where("ATP[c]" == outputMoleculeIDs))[0][0]; ATP = np.sum(deltaMetabolites[:, IdxAtp])
		IdxGtp = (np.where("GTP[c]" == outputMoleculeIDs))[0][0]; GTP = np.sum(deltaMetabolites[:, IdxGtp])
		IdxCtp = (np.where("CTP[c]" == outputMoleculeIDs))[0][0]; CTP = np.sum(deltaMetabolites[:, IdxCtp])
		IdxUtp = (np.where("UTP[c]" == outputMoleculeIDs))[0][0]; UTP = np.sum(deltaMetabolites[:, IdxUtp])
		NtpsProduced = ATP + GTP + CTP + UTP
		ax = plt.subplot(7,2,6)
		plt.plot(time / 60., (deltaMetabolites[:, IdxAtp] + deltaMetabolites[:, IdxGtp] + deltaMetabolites[:, IdxCtp] + deltaMetabolites[:, IdxUtp]) / 1e6)
		plt.ylabel("Metabolism ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("NTPs produced (x$10^{%d}$) = %.2f" % (6, (sum(deltaMetabolites[:, IdxAtp] + deltaMetabolites[:, IdxGtp] + deltaMetabolites[:, IdxCtp] + deltaMetabolites[:, IdxUtp])  / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,9)
		plt.plot(time / 60., FractionActiveEndoRNases * 100)
		plt.ylabel("EndoRN capacity (%)", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,8)
		plt.plot(time / 60., fragmentBasesDigested / 1e6)
		plt.ylabel("Exo-digestion ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("NTPs recycled (x$10^{%d}$) = %.2f" % (6, (fragmentBasesDigested.sum() / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,11)
		plt.plot(time / 60., DiffRelativeFirstOrderDecay)
		plt.ylabel("sum(Residuals)", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		ax = plt.subplot(7,2,10)
		plt.plot(time / 60., (ntpCounts[:, 0] + ntpCounts[:, 1] + ntpCounts[:, 2] + ntpCounts[:, 3]) / 1e6)
		plt.ylabel("Net production ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("NTPs required for cell division (x$10^{%d}$) = %.2f" % (6, ((ntpCounts[0, 0] + ntpCounts[0, 1] + ntpCounts[0, 2] + ntpCounts[0, 3])  / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		# compute active ExoRNase capacity (%)
		ActiveExoRNcapacity = fragmentBasesDigested.astype(float) / (totalexoRnaseCounts * sim_data.constants.KcatExoRNase.asNumber()) * 100

		ax = plt.subplot(7,2,13)
		plt.plot(time / 60., ActiveExoRNcapacity)
		plt.xlabel("Time (min)")
		plt.ylabel("ExoRN capacity (%)", fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)

		# compute instantaneous balance of nTPs
		InstantaneousNTPs = - gtpUsed - countNTPsUSed + (deltaMetabolites[:, IdxAtp] + deltaMetabolites[:, IdxGtp] + deltaMetabolites[:, IdxCtp] + deltaMetabolites[:, IdxUtp]) + fragmentBasesDigested

		ax = plt.subplot(7,2,12)
		plt.plot(time / 60., InstantaneousNTPs / 1e6)
		plt.xlabel("Time (min)")
		plt.ylabel("Balance ($10^{%d}$nt)" % 6, fontsize = 9)
		plt.title("Average instantaneous balance (x$10^{%d}$) = %.4f" % (6, (np.mean(InstantaneousNTPs) / 1e6)), fontsize = 9)
		yloc = plt.MaxNLocator(max_yticks); ax.yaxis.set_major_locator(yloc)


		plt.subplots_adjust(hspace = 0.6, wspace = 0.35)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
	Plot().cli()
