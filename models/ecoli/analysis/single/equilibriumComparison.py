"""
Plot empirical Kd's (from the simulation) and their expected value (from the sim_data)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/24/15
"""

from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

IGNORE_FIRST_PERCENTAGE = 0.1


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load data from KB
		sim_data = cPickle.load(open(simDataFile, "rb"))

		stoichMatrix = sim_data.process.equilibrium.stoichMatrix().astype(np.int64)
		ratesFwd = sim_data.process.equilibrium.ratesFwd
		ratesRev = sim_data.process.equilibrium.ratesRev

		nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		cellDensity = sim_data.constants.cellDensity.asNumber(units.g / units.L)

		moleculeNames = sim_data.process.equilibrium.moleculeNames

		# Load time
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Calculate concentration data
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeIdx = {name: i for i, name in enumerate(bulkMolecules.readAttribute("objectNames"))}
		bulkMoleculeCounts = bulkMolecules.readColumn("counts")

		mass = TableReader(os.path.join(simOutDir, "Mass"))
		cellMass = (units.fg * mass.readColumn("cellMass")).asNumber(units.g)

		cellVolume = cellMass / cellDensity

		fig = plt.figure(figsize = (20, 20))
		rows = 6
		cols = 6
		num_subentries = 3

		for idx in xrange(stoichMatrix.shape[1]):

			grid_loc = idx + 1 + (cols*(num_subentries + 1))*( idx / cols)

			reactantIds = [moleculeNames[x] for x in np.where(stoichMatrix[:, idx] < 0)[0]]
			reactantCoeffs = np.abs(stoichMatrix[stoichMatrix[:, idx] < 0, idx])
			reactantIdxs = np.array([bulkMoleculeIdx[x] for x in reactantIds], dtype = np.int64)
			reactantCounts = bulkMoleculeCounts[:, reactantIdxs]
			reactantConcentrations = reactantCounts / (cellVolume[:, np.newaxis] * nAvogadro)

			productIds = [moleculeNames[x] for x in np.where(stoichMatrix[:, idx] > 0)[0]]
			productCoeffs = np.abs(stoichMatrix[stoichMatrix[:, idx] > 0, idx])
			productIdxs = np.array([bulkMoleculeIdx[x] for x in productIds], dtype = np.int64)
			productCounts = bulkMoleculeCounts[:, productIdxs]
			productConcentrations = productCounts / (cellVolume[:, np.newaxis] * nAvogadro)

			empiricalKd = 10**(np.sum(reactantCoeffs * np.log10(reactantConcentrations), axis = 1) - np.sum(productCoeffs * np.log10(productConcentrations), axis = 1))
			empiricalKd[np.isinf(empiricalKd)] = np.nan
			expectedKd = ratesRev[idx] / ratesFwd[idx]

			ax = plt.subplot(rows*(num_subentries + 2), cols, grid_loc)

			ax.plot(time[1:] / 60., empiricalKd[1:], linewidth=1, label="Empirical K_d")
			ax.plot([time[1] / 60, time[-1] / 60], [expectedKd, expectedKd], linestyle="--")

			bbox = None
			ax.set_title("%s" % productIds, fontsize = 6, bbox = bbox)

			# Sets ticks so that they look pretty
			# ymin = np.nanmin(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
			# ymax = np.nanmax(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
			# if np.any(np.isnan(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])):
			ymin = expectedKd / 2.
			ymax = expectedKd * 2.
			if ymin != ymax:
				ax.set_ylim([ymin, ymax])
				ax.set_yticks([ymin, ymax])
				ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize=6)
			ax.set_xticks([])

			# Plot all reactant concentrations for this reaction
			for reactantIndex in xrange(0,np.amin([reactantConcentrations.shape[1]]+[2])):

				# import ipdb; ipdb.set_trace()

				ax = plt.subplot(rows*(num_subentries + 2), cols, grid_loc + (cols * (reactantIndex + 1)))

				ax.plot(time[1:] / 60., reactantConcentrations[1:,reactantIndex], linewidth=1, label="Reactant concentration", color="g")

				bbox = None
				ax.set_title("%s" % reactantIds[reactantIndex], fontsize = 6, bbox = bbox)

				# Sets ticks so that they look pretty
				# ymin = np.nanmin(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
				# ymax = np.nanmax(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
				# if np.any(np.isnan(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])):
				ymin = np.amin(reactantConcentrations[1:,reactantIndex]*.9)
				ymax = np.amax(reactantConcentrations[1:,reactantIndex]*1.1)
				if ymin != ymax:
					ax.set_ylim([ymin, ymax])
					ax.set_yticks([ymin, ymax])
					ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
				ax.spines['top'].set_visible(False)
				ax.spines['bottom'].set_visible(False)
				ax.xaxis.set_ticks_position('none')
				ax.tick_params(which = 'both', direction = 'out', labelsize=6)
				ax.set_xticks([])

			# Plot all product concentrations for this reaction
			for productIndex in xrange(0,np.amin([productConcentrations.shape[1]]+[2])):
				ax = plt.subplot(rows * (num_subentries + 2), cols, grid_loc + (cols * (productIndex+reactantIndex+2)))

				ax.plot(time[1:] / 60., productConcentrations[1:,productIndex], linewidth=1, label="Product concentration", color="r")

				bbox = None
				ax.set_title("%s" % productIds[productIndex], fontsize = 6, bbox = bbox)

				# Sets ticks so that they look pretty
				# ymin = np.nanmin(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
				# ymax = np.nanmax(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
				# if np.any(np.isnan(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])):
				ymin = np.amin(productConcentrations[1:,productIndex]*.9)
				ymax = np.amax(productConcentrations[1:,productIndex]*1.1)
				if ymin != ymax:
					ax.set_ylim([ymin, ymax])
					ax.set_yticks([ymin, ymax])
					ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
				ax.spines['top'].set_visible(False)
				ax.spines['bottom'].set_visible(False)
				ax.xaxis.set_ticks_position('none')
				ax.tick_params(which = 'both', direction = 'out', labelsize=6)
				ax.set_xticks([])

		# Save
		plt.subplots_adjust(hspace = 1, wspace = 1)

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close("all")

		bulkMolecules.close()


if __name__ == "__main__":
	Plot().cli()
