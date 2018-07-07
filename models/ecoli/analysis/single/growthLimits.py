"""
Plots various effects that may be limiting growth

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/2015
"""

from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))
		nAvogadro = sim_data.constants.nAvogadro
		moleculeIds = sim_data.moleculeGroups.aaIDs
		moleculeIds.append('GTP[c] (translation)')
		moleculeIds.extend(sim_data.moleculeGroups.ntpIds)

		# Load time
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		# Load data
		growthLimitsDataFile = TableReader(os.path.join(simOutDir, "GrowthLimits"))

		# Translation
		gtpPoolSize = growthLimitsDataFile.readColumn("gtpPoolSize")
		gtpRequestSize = growthLimitsDataFile.readColumn("gtpRequestSize")
		gtpAllocated = growthLimitsDataFile.readColumn("gtpAllocated")
		gtpUsed = growthLimitsDataFile.readColumn("gtpUsed")

		aaPoolSize = growthLimitsDataFile.readColumn("aaPoolSize")
		aaRequestSize = growthLimitsDataFile.readColumn("aaRequestSize")
		aaAllocated = growthLimitsDataFile.readColumn("aaAllocated")
		aasUsed = growthLimitsDataFile.readColumn("aasUsed")

		# Transcription
		ntpPoolSize = growthLimitsDataFile.readColumn("ntpPoolSize")
		ntpRequestSize = growthLimitsDataFile.readColumn("ntpRequestSize")
		ntpAllocated = growthLimitsDataFile.readColumn("ntpAllocated")
		ntpUsed = growthLimitsDataFile.readColumn("ntpUsed")


		# Create aggregate
		poolSize = np.hstack((
			aaPoolSize,
			gtpPoolSize.reshape(gtpPoolSize.size,1),
			ntpPoolSize
			))
		requestSize = np.hstack((
			aaRequestSize,
			gtpRequestSize.reshape(gtpPoolSize.size,1),
			ntpRequestSize
			))
		allocated = np.hstack((
			aaAllocated,
			gtpAllocated.reshape(gtpPoolSize.size,1),
			ntpAllocated
			))
		used =  np.hstack((
			aasUsed,
			gtpUsed.reshape(gtpPoolSize.size,1),
			ntpUsed
			))

		growthLimitsDataFile.close()


		# Load ATP usage
		# ATPusageListenerFile = TableReader(os.path.join(simOutDir, "ATPhydrolyzedUsageListener"))
		# atpsHydrolyzed = ATPusageListenerFile.readColumn('atpsHydrolyzed')
		# ATPusageListenerFile.close()


		# bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		# bulkMoleculeIds = bulkMolecules.readAttribute("objectNames")
		# bulkMolecules.close()


		# Plot
		fig = plt.figure(figsize = (11, 11))

		# Plot GTP consumption
		# plt.rc('xtick', labelsize=5)
		# plt.rc('ytick', labelsize=5)
		# plt.plot(time / 60., gtpUsed)
		# plt.plot(time / 60., atpsHydrolyzed)
		# plt.xlabel("Time (min)", fontsize = 5)
		# plt.ylabel("GTP counts", fontsize = 5)
		print "GTP consumption by polypeptide elongation = %d" % np.sum(gtpUsed)
		# print "ATP hydrolyzed for non-growth associated maintenance = %d" % np.sum(atpsHydrolyzed)

		for idx in range(len(moleculeIds)):
			ax = plt.subplot(7,4,idx+1)
			ax.plot(time / 60., poolSize[:,idx], linewidth=2, label="pool size", color='k')
			ax.plot(time / 60., requestSize[:,idx], linewidth=2, label="request size", color='b')
			ax.plot(time / 60., allocated[:,idx], linewidth=2, label="allocated", color='r', linestyle='--')
			ax.plot(time / 60., used[:, idx], linewidth=2, label="used", color="c", linestyle=":")

			# Highlight title if request is greater than pool
			bbox = None
			if (poolSize < requestSize)[:,idx].any() or (allocated < requestSize)[:,idx].any():
				bbox = {'facecolor':'red', 'alpha':0.5, 'pad':10}
			elif (used[:,idx] < allocated[:,idx]).any():
				bbox = {'facecolor':'orange', 'alpha':0.5, 'pad':10}
			ax.set_title(moleculeIds[idx], bbox=bbox)

			# Set ticks so this is all easy to read
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize=10)
			ax.set_xticks([time.min() / 60., time.max() / 60.])
			ax.set_yticks([0., np.max([poolSize[:,idx].max(), requestSize[:,idx].max(), allocated[:,idx].max()])])

		# Create legend
		ax = plt.subplot(7,4,len(moleculeIds) + 1)
		ax.plot(0, 0, linewidth=2, label="pool size", color='k')
		ax.plot(0, 0, linewidth=2, label="request size", color='b')
		ax.plot(0, 0, linewidth=2, label="allocated", color='r', linestyle='--')
		ax.plot(0, 0, linewidth=2, label="used", color="c", linestyle=':')
		ax.legend(loc = 10,prop={'size':10})
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.yaxis.set_ticks_position('none')
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_title("Highlights process under-usage", fontsize=12, bbox={'facecolor':'orange', 'alpha':0.5, 'pad':10})
		ax.set_xlabel("Highlights pool/allocation limit", fontsize=12, bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

		# Save
		plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
