#!/usr/bin/env python
"""
Plot normalized mass fractions of proteome for each protein monomer

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/12/16
"""

from __future__ import division

import argparse
import os

import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils.fitting import normalize
from wholecell.utils import units
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths

import scipy.cluster
import sys

def main(seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	
	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

	allDirs = ap.get_cells()

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))
	nAvogadro = sim_data.constants.nAvogadro

	ids_complexation = sim_data.process.complexation.moleculeNames
	ids_complexation_complexes = [ids_complexation[i] for i in np.where((sim_data.process.complexation.stoichMatrix() == 1).sum(axis = 1))[0]]
	ids_equilibrium = sim_data.process.equilibrium.moleculeNames
	ids_equilibrium_complexes = [ids_equilibrium[i] for i in np.where((sim_data.process.equilibrium.stoichMatrix() == 1).sum(axis = 1))[0]]
	ids_translation = sim_data.process.translation.monomerData["id"].tolist()

	proteomeMWs = sim_data.getter.getMass(ids_translation)

	time = []
	proteomeMassFractions = []
	initialTime = 0

	plt.figure(figsize = (8.5, 11))

	for simDir in allDirs:
		simOutDir = os.path.join(simDir, "simOut")

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		counts = bulkMolecules.readColumn("counts")
		bulkMolecules.close()

		monomerIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_translation], np.int)
		equilibriumIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_equilibrium], np.int)
		equilibriumComplexesIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_equilibrium_complexes], np.int)
		complexationIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_complexation], np.int)
		complexationComplexesIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in ids_complexation_complexes], np.int)

		# Load time
		time = np.append(time, TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime)

		# Load mass data
		massReader = TableReader(os.path.join(simOutDir, "Mass"))
		proteinMass = units.fg * massReader.readColumn("proteinMass")
		massReader.close()

		counts[:, equilibriumIndexes] += np.dot(counts[:, equilibriumComplexesIndexes] * -1, np.matrix.transpose(sim_data.process.equilibrium.stoichMatrixMonomers()))
		counts[:, complexationIndexes] += np.dot(counts[:, complexationComplexesIndexes] * -1, np.matrix.transpose(sim_data.process.complexation.stoichMatrixMonomers()))
		
		# Get mass of proteins in cell
		proteomeMasses = 1. / nAvogadro * counts[:, monomerIndexes] * proteomeMWs
		if proteomeMassFractions == []:
			proteomeMassFractions = proteomeMasses.asNumber(units.fg).T / proteinMass.asNumber(units.fg)
		else:
			proteomeMassFractions = np.concatenate((proteomeMassFractions, proteomeMasses.asNumber(units.fg).T / proteinMass.asNumber(units.fg)), axis = 1)

	# Prevent divide by 0
	proteomeMassFractions += 1e-9

	# Normalize based on mass at early time point and cluster similar traces
	normalizedMassFractions = proteomeMassFractions.T / proteomeMassFractions.T[10, :]
	sys.setrecursionlimit(10000)
	linkage = scipy.cluster.hierarchy.linkage(normalizedMassFractions.T)
	dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot = True)

	nPlots = 16
	nProteins = normalizedMassFractions.shape[1]
	
	for i in range(nPlots):
		ax = plt.subplot(nPlots, 1, i + 1)
		indStart = i * nProteins // nPlots
		indEnd = (i + 1) * nProteins // nPlots
		indexes = dendro["leaves"][indStart:indEnd]

		proteomeSubset = normalizedMassFractions[:, indexes]

		ax.plot(time, proteomeSubset)
		plt.ylabel("Mass Fraction of Proteome", fontsize = 6)

		ymin = np.amin(proteomeSubset * 0.9)
		ymax = np.amax(proteomeSubset * 1.1)
		if ymin != ymax:
			ax.set_ylim([ymin, ymax])
			ax.set_yticks([ymin, ymax])
			ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])

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
