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

	recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
	tfs = sorted(set([x.split("__")[-1] for x in recruitmentColNames if x.split("__")[-1] != "alpha"]))
	argRIndex = [i for i, tf in enumerate(tfs) if tf == "CPLX0-228"][0]

	tfBoundIds = [target + "__CPLX0-228" for target in sim_data.tfToFC["CPLX0-228"].keys()]
	synthProbIds = [target + "[c]" for target in sim_data.tfToFC["CPLX0-228"].keys()]

	plt.figure(figsize = (8.5, 13))
	nRows = 9

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

		# Load data from ribosome data listener
		# ribosomeDataReader = TableReader(os.path.join(simOutDir, "RibosomeData"))
		# nTrpATranslated = ribosomeDataReader.readColumn("numTrpATerminated")
		# ribosomeDataReader.close()

		# Load data from bulk molecules
		bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")

		# Get the concentration of intracellular arg
		argId = ["ARG[c]"]
		argIndex = np.array([bulkMoleculeIds.index(x) for x in argId])
		argCounts = bulkMoleculesReader.readColumn("counts")[:, argIndex].reshape(-1)
		argMols = 1. / nAvogadro * argCounts
		volume = cellMass / cellDensity
		argConcentration = argMols * 1. / volume

		# Get the amount of active argR (that is promoter bound)
		argRActiveId = ["CPLX0-228[c]"]
		argRActiveIndex = np.array([bulkMoleculeIds.index(x) for x in argRActiveId])
		argRActiveCounts = bulkMoleculesReader.readColumn("counts")[:, argRActiveIndex].reshape(-1)

		# Get the amount of inactive argR
		argRInactiveId = ["PC00005[c]"]
		argRInactiveIndex = np.array([bulkMoleculeIds.index(x) for x in argRInactiveId])
		argRInactiveCounts = bulkMoleculesReader.readColumn("counts")[:, argRInactiveIndex].reshape(-1)

		# Get the amount of monomeric argR
		argRMonomerId = ["PD00194[c]"]
		argRMonomerIndex = np.array([bulkMoleculeIds.index(x) for x in argRMonomerId])
		argRMonomerCounts = bulkMoleculesReader.readColumn("counts")[:, argRMonomerIndex].reshape(-1)

		# Get the promoter-bound status for all regulated genes
		tfBoundIndex = np.array([bulkMoleculeIds.index(x) for x in tfBoundIds])
		tfBoundCounts = bulkMoleculesReader.readColumn("counts")[:, tfBoundIndex]

		# Get the amount of monomeric carA
		carAProteinId = ["CARBPSYN-SMALL[c]"]
		carAProteinIndex = np.array([bulkMoleculeIds.index(x) for x in carAProteinId])
		carAProteinCounts = bulkMoleculesReader.readColumn("counts")[:, carAProteinIndex].reshape(-1)

		# Get the amount of complexed carA
		carAComplexId = ["CARBPSYN-CPLX[c]"]
		carAComplexIndex = np.array([bulkMoleculeIds.index(x) for x in carAComplexId])
		carAComplexCounts = bulkMoleculesReader.readColumn("counts")[:, carAComplexIndex].reshape(-1)

		# Get the amount of carA mRNA
		carARnaId = ["EG10134_RNA[c]"]
		carARnaIndex = np.array([bulkMoleculeIds.index(x) for x in carARnaId])
		carARnaCounts = bulkMoleculesReader.readColumn("counts")[:, carARnaIndex].reshape(-1)

		bulkMoleculesReader.close()

		# Compute total counts and concentration of carA in monomeric and complexed form
		# (we know the stoichiometry)
		carAProteinTotalCounts = carAProteinCounts + 2 * carAComplexCounts
		carAProteinTotalMols = 1. / nAvogadro * carAProteinTotalCounts
		carAProteinTotalConcentration = carAProteinTotalMols * 1. / volume

		# Compute concentration of carA mRNA
		carARnaMols = 1. / nAvogadro * carARnaCounts
		carARnaConcentration = carARnaMols * 1. / volume

		# Compute the carA mass in the cell
		carAMw = sim_data.getter.getMass(carAProteinId)
		carAMass = 1. / nAvogadro * carAProteinTotalCounts * carAMw

		# Compute the proteome mass fraction
		proteomeMassFraction = carAMass.asNumber(units.fg) / proteinMass.asNumber(units.fg)

		# Get the synthesis probability for all regulated genes
		rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))

		rnaIds = rnaSynthProbReader.readAttribute("rnaIds")
		synthProbIndex = np.array([rnaIds.index(x) for x in synthProbIds])
		synthProbs = rnaSynthProbReader.readColumn("rnaSynthProb")[:, synthProbIndex]

		argRBound = rnaSynthProbReader.readColumn("nActualBound")[:,argRIndex]
		
		rnaSynthProbReader.close()

		# Calculate total argR - active, inactive, bound and monomeric
		argRTotalCounts = 6 * (argRActiveCounts + argRInactiveCounts + argRBound) + argRMonomerCounts

		# Compute moving averages
		width = 100

		tfBoundCountsMA = np.array([np.convolve(tfBoundCounts[:,i], np.ones(width) / width, mode = "same")
				for i in range(tfBoundCounts.shape[1])]).T
		synthProbsMA = np.array([np.convolve(synthProbs[:,i], np.ones(width) / width, mode = "same")
				for i in range(synthProbs.shape[1])]).T
		
		##############################################################
		ax = plt.subplot(nRows, 1, 1)
		ax.plot(time, argConcentration.asNumber(units.umol / units.L), color = "b")
		plt.ylabel("Internal ARG Conc. [uM]", fontsize = 6)

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
		ax = plt.subplot(nRows, 1, 2)
		ax.semilogy(time, argRActiveCounts, color = "b")
		ax.semilogy(time, argRInactiveCounts, color = "r")
		ax.semilogy(time, argRTotalCounts, color = "g")
		plt.ylabel("ArgR Counts", fontsize = 6)
		plt.legend(["Active (hexamer)", "Inactive (hexamer)", "Total (monomeric)"], fontsize = 6)

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
		ax = plt.subplot(nRows, 1, 3)
		ax.plot(time, tfBoundCountsMA)
		plt.ylabel("ArgR Bound To Promoters\n(Moving Average)", fontsize = 6)

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
		ax = plt.subplot(nRows, 1, 4)
		ax.plot(time, synthProbsMA)
		plt.ylabel("Regulated Gene Synthesis Prob.\n(Moving Average)", fontsize = 6)

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
		ax = plt.subplot(nRows, 1, 5)
		ax.plot(time, carAProteinTotalCounts, color = "b")
		plt.ylabel("CarA Counts", fontsize = 6)

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
		ax = plt.subplot(nRows, 1, 6)
		ax.plot(time, carAProteinTotalConcentration.asNumber(units.umol / units.L), color = "b")
		plt.ylabel("CarA Concentration", fontsize = 6)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.2f" % ymin, "%0.2f" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks([])
		##############################################################

		##############################################################
		ax = plt.subplot(nRows, 1, 7)
		ax.plot(time, carARnaCounts, color = "b")
		plt.ylabel("CarA mRNA Counts", fontsize = 6)

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
		ax = plt.subplot(nRows, 1, 8)
		ax.plot(time, carARnaConcentration.asNumber(units.umol / units.L), color = "b")
		plt.ylabel("CarA mRNA\nConcentration", fontsize = 6)

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
		ax = plt.subplot(nRows, 1, 9)
		ax.plot(time / 3600., proteomeMassFraction, color = "b")
		plt.ylabel("CarA Mass Fraction\nof Proteome", fontsize = 6)

		ymin, ymax = ax.get_ylim()
		ax.set_yticks([ymin, ymax])
		ax.set_yticklabels(["%0.2e" % ymin, "%0.2e" % ymax])
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		# ax.xaxis.set_ticks_position('none')
		ax.tick_params(which = 'both', direction = 'out', labelsize = 6)
		ax.set_xticks(ax.get_xlim())
		##############################################################

	plt.subplots_adjust(hspace = 0.5)

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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])
