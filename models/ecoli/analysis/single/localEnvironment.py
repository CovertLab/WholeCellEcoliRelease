#!/usr/bin/env python
'''
Plots environment nutrient concentrations

@organization: Covert Lab, Department of Bioengineering, Stanford University
'''

from __future__ import absolute_import

import os
import cPickle
import math

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools

from wholecell.io.tablereader import TableReader

from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
from wholecell.analysis.plotting_tools import COLORS_LARGE
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

RANGE_THRESHOLD = 2
MOVING_AVE_WINDOW_SIZE = 200

# set upper end on time (sec). False disables this feature
XLIMIT = False #2000

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, 'simOutDir does not currently exist as a directory'

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load masses
		mass = TableReader(os.path.join(simOutDir, 'Mass'))
		protein = mass.readColumn('proteinMass')
		tRna = mass.readColumn('tRnaMass')
		rRna = mass.readColumn('rRnaMass')
		mRna = mass.readColumn('mRnaMass')
		dna = mass.readColumn('dnaMass')
		smallMolecules = mass.readColumn('smallMoleculeMass')
		total_dry_mass = mass.readColumn('dryMass')

		mass_fractions = np.vstack([
			protein / protein[0],
			rRna / rRna[0],
			tRna / tRna[0],
			mRna / mRna[0],
			dna / dna[0],
			smallMolecules / smallMolecules[0],
			]).T
		massLabels = ['Protein', 'rRNA', 'tRNA', 'mRNA', 'DNA', 'Small Mol.s']

		# Load time
		initialTime = TableReader(
			os.path.join(simOutDir, 'Main')).readAttribute('initialTime')
		time = TableReader(os.path.join(simOutDir, 'Main')).readColumn(
			'time') - initialTime

		# Load environment data
		environment = TableReader(os.path.join(simOutDir, 'Environment'))
		nutrient_names = environment.readAttribute('objectNames')
		nutrient_concentrations = environment.readColumn('nutrientConcentrations')
		nutrient_condition = environment.readColumn('nutrientCondition')
		environment.close()

		# get times of nutrient shifts
		nutrient_shift = list(nutrient_condition)
		nutrient_shift_indices = [0]
		index = 1
		i = 1
		while i < len(nutrient_shift):
			if nutrient_shift[i] == nutrient_shift[i - 1]:
				nutrient_shift.pop(i)
				i -= 1
			else:
				nutrient_shift_indices.append(index)
			index += 1
			i += 1
		nutrient_shift_times = [time[idx] for idx in nutrient_shift_indices]

		# remove blank spaces
		nutrient_shift = [c.replace(' ', '') for c in nutrient_shift]

		# Load flux data
		fbaResults = TableReader(os.path.join(simOutDir, 'FBAResults'))
		reactionIDs = np.array(fbaResults.readAttribute('reactionIDs'))
		externalExchangeFluxes = np.array(fbaResults.readColumn('externalExchangeFluxes'))
		externalExchangeMolecules = np.array(fbaResults.readAttribute('externalExchangeMolecules'))
		importConstraints = np.array(fbaResults.readColumn('importConstraint'))
		importExchanges = np.array(fbaResults.readColumn('importExchange'))
		fbaResults.close()

		# Build a mapping from nutrient_name to color
		idToColor = {}
		for nutrient_name, color in itertools.izip(nutrient_names, itertools.cycle(COLORS_LARGE)):
			idToColor[nutrient_name] = color

		# Build a mapping from reactionID to color
		rxnIdToColor = {}
		for reactionID, color in itertools.izip(reactionIDs, itertools.cycle(COLORS_LARGE)):
			rxnIdToColor[reactionID] = color

		fig = plt.figure(figsize=(30, 30))
		ax1_1 = plt.subplot2grid((7, 2), (0, 0), rowspan=1, colspan=1)
		ax1_2 = plt.subplot2grid((7, 2), (1, 0), rowspan=1, colspan=1)
		ax1_3 = plt.subplot2grid((7, 2), (2, 0), rowspan=4, colspan=1)
		# ax1_4 = plt.subplot2grid((7, 2), (3, 0), rowspan=3, colspan=1)
		ax2_1 = plt.subplot2grid((7, 2), (0, 1), rowspan=1, colspan=1)
		ax2_2 = plt.subplot2grid((7, 2), (2, 1), rowspan=3, colspan=1)

		# total cell mass
		ax1_1.plot(time, total_dry_mass, linewidth=2)
		ax1_1.ticklabel_format(useOffset=False)
		ax1_1.set_xlabel('Time (sec)')
		ax1_1.set_ylabel('Dry Cell Mass (fg)')

		# cell mass fractions
		ax1_2.plot(time, mass_fractions, linewidth=2)
		ax1_2.set_xlabel('Time (sec)')
		ax1_2.set_ylabel('Mass (normalized by t = 0)')
		ax1_2.set_title('Mass Fractions of Biomass Components')
		ax1_2.legend(massLabels, loc='best')

		# plot whether molecule is import exchanged
		ax1_3.imshow(importExchanges.T, aspect='auto', cmap=plt.cm.binary,
			interpolation='nearest', extent=[0,time[-1],len(importConstraints.T)-0.5,-0.5])
		# ax1_3.imshow(importConstraints.T, aspect='auto', cmap=plt.cm.Reds,
		# 	interpolation='nearest', extent=[0,time[-1],len(importConstraints.T)-0.5,-0.5])
		ax1_3.set_title('Import Exchange Molecules (boolean)')
		ax1_3.set_xlabel('Time (sec)')
		ax1_3.set_yticks(np.arange(len(importExchanges.T)))
		ax1_3.set_yticklabels(externalExchangeMolecules, fontsize=8)

		# exchange fluxes
		for idx, (reactionID, externalExchangeFlux) in enumerate(zip(reactionIDs, externalExchangeFluxes.T)):
			runningMeanFlux = np.convolve(externalExchangeFlux, np.ones((MOVING_AVE_WINDOW_SIZE,))/MOVING_AVE_WINDOW_SIZE, mode='valid')
			meanNormFlux = runningMeanFlux / np.mean(runningMeanFlux)
			fluxRange = meanNormFlux.max() - meanNormFlux.min()

			if fluxRange > RANGE_THRESHOLD:
				ax2_1.plot(time, externalExchangeFlux, linewidth=2, label=reactionID, color=rxnIdToColor[reactionID])

		ax2_1.set_yscale('symlog')
		ax2_1.set_xlabel('Time (sec)')
		ax2_1.set_ylabel('symlog Flux {}'.format(FLUX_UNITS.strUnit()))
		ax2_1.set_title('Exchange Fluxes')
		ax2_1.legend(bbox_to_anchor=(0.5, -0.25), loc=9, borderaxespad=0., ncol=1, prop={'size': 10})


		## environment concentrations

		# plot concentrations in environment
		for idx, nutrient_name in enumerate(nutrient_names):
			if (not math.isnan(nutrient_concentrations[0, idx]) and (np.mean(
					nutrient_concentrations[:, idx]) != 0)) and (not math.isinf(nutrient_concentrations[0, idx])):
				ax2_2.plot(time, nutrient_concentrations[:, idx], linewidth=2,
					label=nutrient_name, color=idToColor[nutrient_name])
		ax2_2.set_yscale('symlog',linthreshy=10, linscaley=1)
		ax2_2.set_title('Environment Concentrations')
		ax2_2.set_xlabel('Time (sec)')
		ax2_2.set_ylabel('symlog concentration (mmol/L)')
		ax2_2.legend(bbox_to_anchor=(0.5, -0.25), loc=9, borderaxespad=0.,
			ncol=3, prop={'size': 10})

		if XLIMIT:
			ax1_1.set_xlim([0, XLIMIT])
			ax1_2.set_xlim([0, XLIMIT])
			ax1_3.set_xlim([0, XLIMIT])
			# ax1_4.set_xlim([0, XLIMIT])
			ax2_1.set_xlim([0, XLIMIT])
			ax2_2.set_xlim([0, XLIMIT])

		else:
			ax1_1.set_xlim([0, time[-1]])
			ax1_2.set_xlim([0, time[-1]])
			# ax1_3.set_xlim([0, time[-1]])
			# ax1_4.set_xlim([0, time[-1]])
			ax2_1.set_xlim([0, time[-1]])
			ax2_2.set_xlim([0, time[-1]])


		# plot vertical lines at nutrient shifts
		ymin, ymax = ax1_1.get_ylim()
		for time in nutrient_shift_times:
			ax1_1.plot([time, time], [ymin, ymax], color='k', linestyle='--')

		ymin, ymax = ax1_2.get_ylim()
		for time in nutrient_shift_times:
			ax1_2.plot([time, time], [ymin, ymax], color='k', linestyle='--')

		# ymin, ymax = ax1_3.get_ylim()
		# for time in nutrient_shift_times:
		# 	ax1_3.plot([time, time], [ymin, ymax], color='k', linestyle='--')

		# ymin, ymax = ax1_4.get_ylim()
		# for time in nutrient_shift_times:
		# 	ax1_4.plot([time, time], [ymin, ymax], color='k', linestyle='--')

		ymin, ymax = ax2_1.get_ylim()
		for time in nutrient_shift_times:
			ax2_1.plot([time, time], [ymin, ymax], color='k', linestyle='--')

		ymin, ymax = ax2_2.get_ylim()
		for time in nutrient_shift_times:
			ax2_2.plot([time, time], [ymin, ymax], color='k', linestyle='--')

		# add text for nutrient shifts
		ymin, ymax = ax1_1.get_ylim()
		for idx, (nutrient_condition, time) in enumerate(zip(nutrient_shift, nutrient_shift_times)):
			ax1_1.text(time, ymax+1, nutrient_condition, ha='left', va='bottom', rotation=45)


		plt.subplots_adjust(wspace=0.2, hspace=0.4)

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close('all')



if __name__ == '__main__':
	Plot().cli()