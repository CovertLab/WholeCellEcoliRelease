"""
Plots various simulation components that may be limiting growth

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/2015
"""

from __future__ import absolute_import
from __future__ import division

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import filepath


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		moleculeIds = sim_data.moleculeGroups.aaIDs
		moleculeIds.append('GTP[c] (translation)')
		moleculeIds.extend(sim_data.moleculeGroups.ntpIds)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))
		growth_limits_reader = TableReader(os.path.join(simOutDir, "GrowthLimits"))

		# Load time
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Translation
		gtpPoolSize = growth_limits_reader.readColumn("gtpPoolSize")
		gtpRequestSize = growth_limits_reader.readColumn("gtpRequestSize")
		gtpAllocated = growth_limits_reader.readColumn("gtpAllocated")
		gtpUsed = growth_limits_reader.readColumn("gtpUsed")

		aaPoolSize = growth_limits_reader.readColumn("aaPoolSize")
		aaRequestSize = growth_limits_reader.readColumn("aaRequestSize")
		aaAllocated = growth_limits_reader.readColumn("aaAllocated")
		aasUsed = growth_limits_reader.readColumn("aasUsed")

		# Transcription
		ntpPoolSize = growth_limits_reader.readColumn("ntpPoolSize")
		ntpRequestSize = growth_limits_reader.readColumn("ntpRequestSize")
		ntpAllocated = growth_limits_reader.readColumn("ntpAllocated")
		ntpUsed = growth_limits_reader.readColumn("ntpUsed")

		# Create aggregate
		poolSize = np.hstack((
			aaPoolSize,
			gtpPoolSize.reshape(gtpPoolSize.size, 1),
			ntpPoolSize,
			))
		requestSize = np.hstack((
			aaRequestSize,
			gtpRequestSize.reshape(gtpPoolSize.size, 1),
			ntpRequestSize,
			)).astype(np.int64)
		allocated = np.hstack((
			aaAllocated,
			gtpAllocated.reshape(gtpPoolSize.size, 1),
			ntpAllocated,
			))
		used =  np.hstack((
			aasUsed,
			gtpUsed.reshape(gtpPoolSize.size, 1),
			ntpUsed,
			))

		# Plot
		fig = plt.figure(figsize = (11, 11))

		for idx in range(len(moleculeIds)):
			ax = plt.subplot(7, 4, idx+1)
			ax.plot(time / 60., poolSize[:,idx], linewidth=2, label="pool size", color='k')
			ax.plot(time / 60., requestSize[:,idx], linewidth=2, label="request size", color='b')
			ax.plot(time / 60., allocated[:,idx], linewidth=2, label="allocated", color='r', linestyle='--')
			ax.plot(time / 60., used[:, idx], linewidth=2, label="used", color="c", linestyle=":")

			# Highlight title if request is greater than pool
			bbox = None
			if (poolSize[:, idx] < requestSize[:, idx]).any() or (allocated[:, idx] < requestSize[:, idx]).any():
				bbox = {'facecolor':'red', 'alpha':0.5, 'pad':10}
			elif (used[:, idx] < allocated[:, idx]).any():
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
		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
