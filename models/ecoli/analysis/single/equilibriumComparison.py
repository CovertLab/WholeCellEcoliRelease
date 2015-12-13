#!/usr/bin/env python
"""
Plot empirical Kd's (from the simulation) and their expected value (from the sim_data)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/24/15
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

IGNORE_FIRST_PERCENTAGE = 0.1

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))
	validation_data = cPickle.load(open(validationDataFile, "rb"))

	stoichMatrix = sim_data.process.equilibrium.stoichMatrix().astype(np.int64)
	ratesFwd = sim_data.process.equilibrium.ratesFwd
	ratesRev = sim_data.process.equilibrium.ratesRev

	nAvogadro = sim_data.constants.nAvogadro.asNumber(1 / units.mol)
	cellDensity = sim_data.constants.cellDensity.asNumber(units.g / units.L)


	moleculeNames = sim_data.process.equilibrium.moleculeNames


	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Calculate concentration data
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	bulkMoleculeIds = bulkMolecules.readAttribute("objectNames")

	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = (units.fg * mass.readColumn("cellMass")).asNumber(units.g)

	cellVolume = cellMass / cellDensity

	fig = plt.figure(figsize = (20, 20))
	rows = 7
	cols = 3
	num_subentries = 3

	for idx in xrange(stoichMatrix.shape[1]):

		grid_loc = idx + 1 + (cols*(num_subentries + 1))*( idx / cols)

		reactantIds = [moleculeNames[x] for x in np.where(stoichMatrix[:, idx] < 0)[0]]
		reactantIdxs = np.array([bulkMoleculeIds.index(x) for x in reactantIds], dtype = np.int64)
		reactantCounts = bulkMolecules.readColumn("counts")[:, reactantIdxs]
		reactantConcentrations = reactantCounts / (cellVolume[:, np.newaxis] * nAvogadro)

		productIds = [moleculeNames[x] for x in np.where(stoichMatrix[:, idx] > 0)[0]]
		productIdxs = np.array([bulkMoleculeIds.index(x) for x in productIds], dtype = np.int64)
		productCounts = bulkMolecules.readColumn("counts")[:, productIdxs]
		productConcentrations = productCounts / (cellVolume[:, np.newaxis] * nAvogadro)

		empiricalKd = reactantConcentrations.prod(axis = 1) / productConcentrations.prod(axis = 1)
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

			ax = plt.subplot(rows*(num_subentries + 2), cols, grid_loc+((cols)*(reactantIndex+1)))

			ax.plot(time[1:] / 60., reactantConcentrations[1:,reactantIndex], linewidth=1, label="Reactant concentration", color="g")

			bbox = None
			ax.set_title("%s" % reactantIds[reactantIndex], fontsize = 6, bbox = bbox)

			# Sets ticks so that they look pretty
			# ymin = np.nanmin(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
			# ymax = np.nanmax(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
			# if np.any(np.isnan(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])):
			ymin = np.amin(reactantConcentrations[1:,reactantIndex]*.9)
			ymax = np.amax(reactantConcentrations[1:,reactantIndex]*1.1)
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
			ax = plt.subplot(rows*(num_subentries + 2), cols, grid_loc+((cols)*(productIndex+reactantIndex+2)))

			ax.plot(time[1:] / 60., productConcentrations[1:,productIndex], linewidth=1, label="Product concentration", color="r")

			bbox = None
			ax.set_title("%s" % productIds[productIndex], fontsize = 6, bbox = bbox)

			# Sets ticks so that they look pretty
			# ymin = np.nanmin(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
			# ymax = np.nanmax(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])
			# if np.any(np.isnan(empiricalKd[empiricalKd.shape[0] * IGNORE_FIRST_PERCENTAGE:])):
			ymin = np.amin(productConcentrations[1:,productIndex]*.9)
			ymax = np.amax(productConcentrations[1:,productIndex]*1.1)
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
			ax.spines['top'].set_visible(False)
			ax.spines['bottom'].set_visible(False)
			ax.xaxis.set_ticks_position('none')
			ax.tick_params(which = 'both', direction = 'out', labelsize=6)
			ax.set_xticks([])


	# Create legend
	# ax = plt.subplot(rows, cols, stoichMatrix.shape[1] + 2)
	# ax.plot(0, 0, linewidth=2, label="count", color='k')
	# ax.legend(loc = 10,prop={'size':10})
	# ax.spines['top'].set_visible(False)
	# ax.spines['bottom'].set_visible(False)
	# ax.spines['left'].set_visible(False)
	# ax.spines['right'].set_visible(False)
	# ax.xaxis.set_ticks_position('none')
	# ax.yaxis.set_ticks_position('none')
	# ax.set_xticks([])
	# ax.set_yticks([])
	# ax.set_title("Highlights low empiricalKd", fontsize=12, bbox={'facecolor':'red', 'alpha':0.5, 'pad':1})

	# Save
	plt.subplots_adjust(hspace = 1, wspace = 1)

	from wholecell.analysis.analysis_tools import exportFigure

	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

	plt.close("all")

	bulkMolecules.close()

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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])
