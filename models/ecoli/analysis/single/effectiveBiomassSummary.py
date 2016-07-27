#!/usr/bin/env python
"""
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/25/2016
"""

from __future__ import division

import argparse
import os
import json
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS
FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS

BURN_IN_PERIOD = 100

MEAN_DIFFERENCE_TOLERANCE = 5
STD_DIFFERENCE_TOLERANCE = 5

MEAN_AND_STD = 'purple'
MEAN_ONLY = 'red'
STD_ONLY = 'blue'
EXPECTED_VALUE_COLOR = 'black'
EXPECTED_LOG_VALUE_COLOR = 'green'

EXCEPTION_ON_FLUX_DIFF = False

LOWER_CONC_RELEVANCE_BOUND = 1e-10

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	outputFluxes = fbaResults.readColumn("outputFluxes")
	outputMoleculeIDs = np.array(fbaResults.readAttribute("outputMoleculeIDs"))
	fbaResults.close()

	mean_biomass = np.mean(outputFluxes[BURN_IN_PERIOD:], axis=0)
	mean_log10_biomass = np.mean(np.log10(np.abs(outputFluxes[BURN_IN_PERIOD:])), axis=0)
	std_biomass = np.std(outputFluxes[BURN_IN_PERIOD:], axis=0)

	mean_log10_biomass[mean_log10_biomass == -np.inf] = -25

	sim_data = cPickle.load(open(simDataFile, "rb"))

	outputFluxesMod = outputFluxes
	outputMoleculeIDsMod = outputMoleculeIDs

	fig = plt.figure(figsize = (30, 15))

	data, logData, meanNormData, logMeanNormData = [], [], [], []
	for molIdx, (molID, timecourse) in enumerate(zip(outputMoleculeIDsMod, outputFluxesMod.T)):
		data.append(timecourse[BURN_IN_PERIOD:])
		logData.append(np.log10(timecourse)[BURN_IN_PERIOD:])
		meanNormData.append(timecourse[BURN_IN_PERIOD:] / np.mean(timecourse[BURN_IN_PERIOD:]))
		logMeanNormData.append(np.log10(timecourse[BURN_IN_PERIOD:] / np.mean(timecourse[BURN_IN_PERIOD:])))

	# Read previous biomass fluxes from file
	previousBiomassMeans = {key:value.asNumber(COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) for key, value in sim_data.process.metabolism.previousBiomassMeans.iteritems()}
	previousBiomassLog10Means = {key:value.asNumber(COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) for key, value in sim_data.process.metabolism.previousBiomassLog10Means.iteritems()}
	previousBiomassStds = {key:value.asNumber(COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS) for key, value in sim_data.process.metabolism.previousBiomassStds.iteritems()}

	# Find output fluxes differing by more than a given factor from the values predicted in reconstruction
	differingMeans, differingStds = set(), set()
	for molID, simulationMean, simulationStd in zip(outputMoleculeIDs, mean_biomass, std_biomass):
		if molID not in previousBiomassMeans or molID not in previousBiomassStds:
			continue
		previousMean = previousBiomassMeans[molID]
		previousStd = previousBiomassStds[molID]
		if np.abs(simulationMean - previousMean) / previousMean > (MEAN_DIFFERENCE_TOLERANCE - 1):
			if simulationMean - previousMean > LOWER_CONC_RELEVANCE_BOUND:
				differingMeans.add(molID)
		if np.abs(simulationStd - previousStd) / previousStd > (STD_DIFFERENCE_TOLERANCE - 1):
			if simulationStd - previousStd > LOWER_CONC_RELEVANCE_BOUND:
				differingStds.add(molID)

	plt.suptitle("Fluxes for Output Molecule Reactions")

	xTicksList = []
	plt.subplot(5,1,1)
	plt.boxplot(data)
	xTicks = plt.xticks(range(1,len(outputMoleculeIDsMod)+1), [x for x in outputMoleculeIDsMod], rotation='vertical', fontsize='x-small')
	plt.ylabel('Output Molecule Flux')
	xTicksList.append(xTicks)

	plt.subplot(5,1,2)
	plt.boxplot(logData)
	xTicks = plt.xticks(range(1,len(outputMoleculeIDsMod)+1), [x for x in outputMoleculeIDsMod], rotation='vertical', fontsize='x-small')
	plt.ylabel('Log10')
	xTicksList.append(xTicks)

	plt.subplot(5,1,3)
	plt.boxplot(meanNormData)
	xTicks = plt.xticks(range(1,len(outputMoleculeIDsMod)+1), [x for x in outputMoleculeIDsMod], rotation='vertical', fontsize='x-small')
	plt.ylabel('Mean Normalized')
	xTicksList.append(xTicks)

	plt.subplot(5,1,4)
	plt.boxplot(logMeanNormData)
	xTicks = plt.xticks(range(1,len(outputMoleculeIDsMod)+1), [x for x in outputMoleculeIDsMod], rotation='vertical', fontsize='x-small')
	plt.ylabel('Log10 Mean Normalized')
	xTicksList.append(xTicks)

	for idx, outputMoleculeID in enumerate(outputMoleculeIDsMod):
		if outputMoleculeID in differingMeans and outputMoleculeID in differingStds:
			statusColor = MEAN_AND_STD
		elif outputMoleculeID in differingMeans:
			statusColor = MEAN_ONLY
		elif outputMoleculeID in differingStds:
			statusColor = STD_ONLY
		else:
			continue
		
		for xTicks in xTicksList:
			xTicks[0][idx].label.set_color(statusColor)

	textArea = plt.subplot(5,1,5)
	textArea.axis('off')
	plt.ylim([0,10])
	plt.text(0,1,'{} text indicates mean is at least {} fold different from reconstruction flat file.'.format(MEAN_ONLY, MEAN_DIFFERENCE_TOLERANCE),color=MEAN_ONLY)
	plt.text(0,0,'{} text indicates std is at least {} fold different from reconstruction flat file.'.format(STD_ONLY, STD_DIFFERENCE_TOLERANCE),color=STD_ONLY)
	plt.text(0,-1,'{} text indicates both mean and std are at least {} and {} fold different from reconstruction flat file, respectively.'.format(MEAN_AND_STD, MEAN_DIFFERENCE_TOLERANCE, STD_DIFFERENCE_TOLERANCE),color=MEAN_AND_STD)
	plt.text(0,-3,'{} horizontal line indicates the value expected from the reconstruction flat file.'.format(EXPECTED_VALUE_COLOR),color=EXPECTED_VALUE_COLOR)
	plt.text(0,-4,'{} horizontal line indicates the mean of the log10 value expected from the reconstruction flat file.'.format(EXPECTED_LOG_VALUE_COLOR),color=EXPECTED_LOG_VALUE_COLOR)
	plt.text(0,-5,'red horizontal line indicates the simulation mean flux.',color='red')


	outputMoleculeIDsModList = list(outputMoleculeIDsMod)
	for outputMoleculeID, mean in previousBiomassMeans.iteritems():
		idx = outputMoleculeIDsModList.index(outputMoleculeID)
		plt.subplot(5,1,1)
		plt.hlines(y=mean,xmin=idx+.5,xmax=idx+1.5,color=EXPECTED_VALUE_COLOR)
		plt.subplot(5,1,2)
		plt.hlines(y=np.log10(np.abs(mean)),xmin=idx+.5,xmax=idx+1.5,color=EXPECTED_VALUE_COLOR)
		plt.hlines(y=previousBiomassLog10Means[outputMoleculeID],xmin=idx+.5,xmax=idx+1.5,color=EXPECTED_LOG_VALUE_COLOR)


	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")


	# Write current biomass fluxes to file
	with open(os.path.join(plotOutDir, plotOutFileName + '.tsv'), 'w') as output:
		output.write(
			"\"molecule id\"\t\"mean flux (units.{counts} / units.{volume} / units.{time})\"\t\"mean log10 flux (units.{counts} / units.{volume} / units.{time})\"\t\"standard deviation (units.{counts} / units.{volume} / units.{time})\"\n".format(
				counts=COUNTS_UNITS.strUnit().strip('[]'),
				volume=VOLUME_UNITS.strUnit().strip('[]'),
				time=TIME_UNITS.strUnit().strip('[]')
			)
		)
		for molID, meanFlux, meanLogFlux, stdFlux in zip(outputMoleculeIDs, mean_biomass, mean_log10_biomass, std_biomass):
			output.write("\t".join(['"'+unicode.encode(molID)+'"', str(meanFlux), str(meanLogFlux), str(stdFlux)]) + "\n")

	# Throw an exception if this cell's biomass fluxes do not match the flat file
	# (indicates that the flat file needs to be changed/updated)
	if (len(differingMeans) > 0 or len(differingStds) > 0):
		if EXCEPTION_ON_FLUX_DIFF:
			raise Exception(
				"Biomass fluxes predicted from the flat file in reconstruction do not match simulation fluxes. {} means differ by a factor of {} or more: {}\n\n, {} std differ by a factor of {} or more: {}".format(
					len(differingMeans), MEAN_DIFFERENCE_TOLERANCE, differingMeans, len(differingStds), STD_DIFFERENCE_TOLERANCE, differingStds)
				)
		else:
			print "Biomass fluxes predicted from the flat file in reconstruction do not match simulation fluxes. {} means differ by a factor of {} or more: {}\n\n, {} std differ by a factor of {} or more: {}".format(
					len(differingMeans), MEAN_DIFFERENCE_TOLERANCE, differingMeans, len(differingStds), STD_DIFFERENCE_TOLERANCE, differingStds)


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