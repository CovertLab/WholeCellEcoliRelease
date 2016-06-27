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

MAX_LEN = 10
SAMPLE_EVERY = 1
BURN_IN_PERIOD = 100

NUMERICAL_ZERO = 1e-9

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	outputFluxes = fbaResults.readColumn("outputFluxes")
	outputMoleculeIDs = np.array(fbaResults.readAttribute("outputMoleculeIDs"))
	fbaResults.close()

	mean_biomass = np.mean(outputFluxes, axis=0)
	std_biomass = np.std(outputFluxes, axis=0)

	sim_data = cPickle.load(open(simDataFile, "rb"))

	# Read previous biomass fluxes from file
	observedBiomassMeans = sim_data.process.metabolism.observedBiomassMeans
	observedBiomassStds = sim_data.process.metabolism.observedBiomassStds

	# Throw an exception if this cell's biomass fluxes do not match the flat file
	# (indicates that the flat file needs to be changed/updated)
	differingMeans, differingStds = set(), set()
	for molID, simulationMean, simulationStd in zip(outputMoleculeIDs, mean_biomass, std_biomass):
		if molID not in observedBiomassMeans or molID not in observedBiomassStds:
			continue
		observedMean = observedBiomassMeans[molID]
		observedStd = observedBiomassStds[molID]
		if np.abs(simulationMean - observedMean) > NUMERICAL_ZERO:
			differingMeans.add(molID)
		if np.abs(simulationStd - observedStd) > NUMERICAL_ZERO:
			differingStds.add(molID)
	if len(differingMeans) > 0 or len(differingStds) > 0:
		raise Exception(
			"Biomass fluxes predicted from the flat file in reconstruction do not match simulation fluxes. Differing means: {}, Differing std: {}".format(
				differingMeans, differingStds)
			)

	# Write current biomass fluxes to file
	with open(os.path.join(plotOutDir, plotOutFileName + '.tsv'), 'w') as output:
		output.write("\"molecule id\"\t\"mean flux\"\t\"standard deviation\"\n")
		for molID, meanFlux, stdFlux in zip(outputMoleculeIDs, mean_biomass, std_biomass):
			output.write("\t".join(['"'+unicode.encode(molID)+'"', str(meanFlux), str(stdFlux)]) + "\n")

	# maxFluxIdx = np.argmax(mean_biomass)

	# outputFluxesMod = np.concatenate((outputFluxes[:,:maxFluxIdx], outputFluxes[:,maxFluxIdx+1:]),axis=1)
	# outputMoleculeIDsMod = np.concatenate((outputMoleculeIDs[:maxFluxIdx], outputMoleculeIDs[maxFluxIdx+1:]),axis=1)

	outputFluxesMod = outputFluxes
	outputMoleculeIDsMod = outputMoleculeIDs

	fig = plt.figure(figsize = (30, 15))

	data, logData, meanNormData = [], [], []
	for molIdx, (molID, timecourse) in enumerate(zip(outputMoleculeIDsMod, outputFluxesMod.T)):
		data.append(timecourse[BURN_IN_PERIOD::SAMPLE_EVERY])
		logData.append(np.log10(timecourse)[BURN_IN_PERIOD::SAMPLE_EVERY])
		meanNormData.append(timecourse / np.mean(timecourse))

	plt.subplot(4,1,1)
	plt.boxplot(data)
	plt.xticks(range(1,len(outputMoleculeIDsMod)), [x[:-3][:MAX_LEN] for x in outputMoleculeIDsMod], rotation='vertical', fontsize='x-small')
	plt.ylabel('Output Molecule Flux')

	plt.subplot(4,1,2)
	plt.boxplot(logData)
	plt.xticks(range(1,len(outputMoleculeIDsMod)), [x[:-3] for x in outputMoleculeIDsMod], rotation='vertical', fontsize='x-small')
	plt.ylabel('Log10 Output Molecule Flux')

	plt.subplot(4,1,3)
	plt.boxplot(meanNormData)
	plt.xticks(range(1,len(outputMoleculeIDsMod)), [x[:-3] for x in outputMoleculeIDsMod], rotation='vertical', fontsize='x-small')
	plt.ylabel('Mean Normalized Output Molecule Flux')

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
	plt.close("all")

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
