#!/usr/bin/env python
"""
Plot trp regulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/2016
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
	cellDensity = sim_data.constants.cellDensity

	plt.figure(figsize = (8.5, 11))

	for simDir in allDirs:
		simOutDir = os.path.join(simDir, "simOut")
		# Load time
		initialTime = 0#TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		# Load mass data
		# Total cell mass is needed to compute concentrations (since we have cell density)
		# Protein mass is needed to compute the mass fraction of the proteome that is trpA
		massReader = TableReader(os.path.join(simOutDir, "Mass"))

		cellMass = units.fg * massReader.readColumn("cellMass")
		proteinMass = units.fg * massReader.readColumn("proteinMass")


		massReader.close()

		# Load data from bulk molecules
		bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")

		# Get the concentration of intracellular trp
		trpId = ["TRP[c]"]
		trpIndex = np.array([bulkMoleculeIds.index(x) for x in trpId])
		trpCounts = bulkMoleculesReader.readColumn("counts")[:, trpIndex].reshape(-1)
		trpMols = 1. / nAvogadro * trpCounts
		volume = cellMass / cellDensity
		trpConcentration = trpMols * 1. / volume

		# Get the amount of active trpR (that isn't promoter bound)
		trpRActiveId = ["CPLX-125[c]"]
		trpRActiveIndex = np.array([bulkMoleculeIds.index(x) for x in trpRActiveId])
		trpRActiveCounts = bulkMoleculesReader.readColumn("counts")[:, trpRActiveIndex].reshape(-1)

		# Get the promoter-bound status of the trpA gene
		trpATfBoundId = ["EG11024_RNA__CPLX-125"]
		trpATfBoundIndex = np.array([bulkMoleculeIds.index(x) for x in trpATfBoundId])
		trpATfBoundCounts = bulkMoleculesReader.readColumn("counts")[:, trpATfBoundIndex].reshape(-1)

		# Get the amount of monomeric trpA
		trpAProteinId = ["TRYPSYN-APROTEIN[c]"]
		trpAProteinIndex = np.array([bulkMoleculeIds.index(x) for x in trpAProteinId])
		trpAProteinCounts = bulkMoleculesReader.readColumn("counts")[:, trpAProteinIndex].reshape(-1)

		# Get the amount of complexed trpA
		trpABComplexId = ["TRYPSYN[c]"]
		trpABComplexIndex = np.array([bulkMoleculeIds.index(x) for x in trpABComplexId])
		trpABComplexCounts = bulkMoleculesReader.readColumn("counts")[:, trpABComplexIndex].reshape(-1)

		bulkMoleculesReader.close()

		# Compute total counts of trpA in monomeric and complexed form
		# (we know the stoichiometry)
		trpAProteinTotalCounts = trpAProteinCounts + 2 * trpABComplexCounts

		# Compute the trpA mass in the cell
		trpAMw = sim_data.getter.getMass(trpAProteinId)
		trpAMass = 1. / nAvogadro * trpAProteinTotalCounts * trpAMw


		# Compute the proteome mass fraction
		proteomeMassFraction = trpAMass.asNumber(units.fg) / proteinMass.asNumber(units.fg)

		# Get the trpA synthesis probability
		rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))

		rnaIds = rnaSynthProbReader.readAttribute("rnaIds")
		trpASynthProbId = ["EG11024_RNA[c]"]
		trpASynthProbIndex = np.array([rnaIds.index(x) for x in trpASynthProbId])
		trpASynthProb = rnaSynthProbReader.readColumn("rnaSynthProb")[:, trpASynthProbIndex].reshape(-1)
		
		rnaSynthProbReader.close()

		# Compute moving averages
		width = 100

		trpATfBoundCountsMA = np.convolve(trpATfBoundCounts, np.ones(width) / width, mode = "same")
		trpASynthProbMA = np.convolve(trpASynthProb, np.ones(width) / width, mode = "same")
		
		##############################################################
		ax = plt.subplot(6, 1, 1)
		ax.plot(time, trpConcentration.asNumber(units.umol / units.L), color = "b")
		plt.ylabel("Internal TRP Conc. [uM]", fontsize = 6)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(6, 1, 2)
		ax.plot(time, trpRActiveCounts, color = "b")
		plt.ylabel("Active TrpR Counts", fontsize = 6)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize=6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(6, 1, 3)
		ax.plot(time, trpATfBoundCounts, color = "b")
		ax.plot(time, trpATfBoundCountsMA, color = "g")
		plt.ylabel("TrpR Bound To trpA Promoter", fontsize = 6)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(6, 1, 4)
		ax.plot(time, trpASynthProb, color = "b")
		ax.plot(time, trpASynthProbMA, color = "g")
		plt.ylabel("trpA Synthesis Prob.", fontsize = 6)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(6, 1, 5)
		ax.plot(time, trpAProteinTotalCounts, color = "b")
		plt.ylabel("TrpA Counts", fontsize = 6)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.0f" % ymin, "%0.0f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################


		##############################################################
		ax = plt.subplot(6, 1, 6)
		ax.plot(time / 3600., proteomeMassFraction, color = "b")
		plt.ylabel("TrpA Mass Fraction of Proteome", fontsize = 6)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		# ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks(ax.get_xlim())
		##############################################################

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
