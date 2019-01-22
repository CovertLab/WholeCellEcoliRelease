'''
Shows shadow prices and reduced costs for the FBA objective
Plots all reduced cost and shadow prices on separate subplots
Identifies outliers for each and plots those on separate subplots

@date: Created 2/6/2018
@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
'''

from __future__ import absolute_import
from __future__ import division

import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader

from wholecell.analysis.plotting_tools import COLORS_SMALL
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Read optimization problem results
		fbaResults = TableReader(os.path.join(simOutDir, 'FBAResults'))
		metaboliteNames = np.array(fbaResults.readAttribute('metaboliteNames'))
		reactionNames = np.array(fbaResults.readAttribute('reactionIDs'))
		kineticReactionNames = np.array(fbaResults.readAttribute('kineticTargetFluxNames'))
		shadowPrices = fbaResults.readColumn('shadowPrices')
		reducedCosts = fbaResults.readColumn('reducedCosts')
		fbaResults.close()

		# Read time info from the listener
		mainReader = TableReader(os.path.join(simOutDir, 'Main'))
		time = mainReader.readColumn('time') - mainReader.readAttribute('initialTime')
		mainReader.close()

		# Set properties for and create figure
		# Burn in required because some large spikes/instability in early time points
		burnIn = 20  # timesteps
		nStd = 1  # number of standard deviations away from mean to filter
		plt.figure(figsize = (8.5, 11))

		# Plot all shadow prices
		ax = self.subplot(3, 2, 1)
		ax.set_prop_cycle('color', COLORS_SMALL)
		ax.plot(time[burnIn:], shadowPrices[burnIn:, :])
		ax.set_xlabel('Time (s)', fontsize=8)
		ax.set_ylabel('Shadow Price\n(Metabolites)', fontsize=8)
		ax.tick_params(axis='both', which='major', labelsize=6)

		# Plot all reduced costs
		ax = self.subplot(3, 2, 3)
		ax.set_prop_cycle('color', COLORS_SMALL)
		ax.plot(time[burnIn:], reducedCosts[burnIn:, :])
		ax.set_xlabel('Time (s)', fontsize=8)
		ax.set_ylabel('Reduced Costs\n(Reactions)', fontsize=8)
		ax.tick_params(axis='both', which='major', labelsize=6)

		# Find outliers for shadow prices and plot with legend to identify
		means = np.mean(shadowPrices, axis=0)
		mean = np.mean(means)
		std = np.std(means)
		outliers = np.unique(np.where((shadowPrices[burnIn:, :] > mean + nStd*std) | (shadowPrices[burnIn:, :] < mean - nStd*std))[1])
		idx = outliers[np.argsort(means[outliers])][::-1]

		ax = self.subplot(3, 2, 2)
		ax.set_prop_cycle('color', COLORS_SMALL)
		ax.set_xlabel('Time (s)', fontsize=8)
		ax.set_ylabel('Shadow Price\n(Metabolites)', fontsize=8)
		if len(idx):
			ax.plot(time[burnIn:], shadowPrices[burnIn:, idx])
			ax.legend(metaboliteNames[idx], fontsize=6)
			ax.tick_params(axis='both', which='major', labelsize=6)

		# Find outliers for reduced costs and plot with legend to identify
		means = np.mean(reducedCosts, axis=0)
		mean = np.mean(means)
		std = np.std(means)
		outliers = np.unique(np.where((reducedCosts[burnIn:, :] > mean + nStd*std) | (reducedCosts[burnIn:, :] < mean - nStd*std))[1])
		idx = outliers[np.argsort(means[outliers])][::-1]

		ax = self.subplot(3, 2, 4)
		ax.set_prop_cycle('color', COLORS_SMALL)
		ax.set_xlabel('Time (s)', fontsize=8)
		ax.set_ylabel('Reduced Costs\n(Reactions)', fontsize=8)
		if len(idx):
			shortenedNames = [x[:20] for x in reactionNames[idx]] # prevent long names from covering both subplots
			ax.plot(time[burnIn:], reducedCosts[burnIn:, idx])
			ax.legend(shortenedNames, fontsize=6)
			ax.tick_params(axis='both', which='major', labelsize=6)

		# Find outliers for reduced costs and plot with legend to identify
		means = np.mean(reducedCosts, axis=0)
		mean = np.mean(means)
		std = np.std(means)
		outliers = np.unique(np.where((reducedCosts[burnIn:, :] > mean + nStd*std) | (reducedCosts[burnIn:, :] < mean - nStd*std))[1])
		idx = outliers[np.argsort(means[outliers])][::-1]

		ax = self.subplot(3, 2, 4)
		ax.set_prop_cycle('color', COLORS_SMALL)
		ax.set_xlabel('Time (s)', fontsize=8)
		ax.set_ylabel('Reduced Costs\n(Reactions)', fontsize=8)
		if len(idx):
			shortenedNames = [x[:20] for x in reactionNames[idx]] # prevent long names from covering both subplots
			ax.plot(time[burnIn:], reducedCosts[burnIn:, idx])
			ax.legend(shortenedNames, fontsize=6)
			ax.tick_params(axis='both', which='major', labelsize=6)

		# Check if there are kinetic reactions (FBA kinetics enabled)
		if len(kineticReactionNames):
			kineticsIdx = [np.where(reactionNames == rxn)[0][0] for rxn in kineticReactionNames]
			kineticReactionNames = np.array(reactionNames)[kineticsIdx]
			kineticsReducedCosts = reducedCosts[:, kineticsIdx]

			# Plot all reduced costs for kinetically constrained reactions
			ax = self.subplot(3, 2, 5)
			ax.set_prop_cycle('color', COLORS_SMALL)
			ax.plot(time[burnIn:], kineticsReducedCosts[burnIn:, :])
			ax.set_xlabel('Time (s)', fontsize=8)
			ax.set_ylabel('Reduced Costs\n(Kinetically Constrained Reactions)', fontsize=8)
			ax.tick_params(axis='both', which='major', labelsize=6)

			# Find outliers for reduced costs and plot with legend to identify
			means = np.mean(kineticsReducedCosts, axis=0)
			mean = np.mean(means)
			std = np.std(means)
			outliers = np.unique(np.where((kineticsReducedCosts[burnIn:, :] > mean + nStd*std) | (kineticsReducedCosts[burnIn:, :] < mean - nStd*std))[1])
			idx = outliers[np.argsort(means[outliers])][::-1]

			ax = self.subplot(3, 2, 6)
			ax.set_prop_cycle('color', COLORS_SMALL)
			ax.set_xlabel('Time (s)', fontsize=8)
			ax.set_ylabel('Reduced Costs\n(Kinetically Constrained Reactions)', fontsize=8)
			if len(idx):
				shortenedNames = [x[:20] for x in kineticReactionNames[idx]]  # prevent long names from covering both subplots
				ax.plot(time[burnIn:], kineticsReducedCosts[burnIn:, idx])
				ax.legend(shortenedNames, fontsize=6)
				ax.tick_params(axis='both', which='major', labelsize=6)

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
