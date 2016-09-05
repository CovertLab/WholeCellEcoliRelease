#!/usr/bin/env python
"""
Plot tyr regulation

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
	tyrRIndex = [i for i, tf in enumerate(tfs) if tf == "MONOMER0-162"][0]

	plt.figure(figsize = (8.5, 11))

	for simDir in allDirs:
		simOutDir = os.path.join(simDir, "simOut")
		# Load time
		initialTime = 0
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		# Load mass data
		# Total cell mass is needed to compute concentrations (since we have cell density)
		# Protein mass is needed to compute the mass fraction of the proteome that is tyrA
		massReader = TableReader(os.path.join(simOutDir, "Mass"))

		cellMass = units.fg * massReader.readColumn("cellMass")
		proteinMass = units.fg * massReader.readColumn("proteinMass")


		massReader.close()

		# Load data from bulk molecules
		bulkMoleculesReader = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulkMoleculeIds = bulkMoleculesReader.readAttribute("objectNames")
		# Get the concentration of intracellular phe
		pheId = ["PHE[c]"]
		pheIndex = np.array([bulkMoleculeIds.index(x) for x in pheId])
		pheCounts = bulkMoleculesReader.readColumn("counts")[:, pheIndex].reshape(-1)
		pheMols = 1. / nAvogadro * pheCounts
		volume = cellMass / cellDensity
		pheConcentration = pheMols * 1. / volume

		# Get the amount of active tyrR (that isn't promoter bound)
		tyrRActiveId = ["MONOMER0-162[c]"]
		tyrRActiveIndex = np.array([bulkMoleculeIds.index(x) for x in tyrRActiveId])
		tyrRActiveCounts = bulkMoleculesReader.readColumn("counts")[:, tyrRActiveIndex].reshape(-1)

		# Get the amount of inactive tyrR
		tyrRInactiveId = ["PD00413[c]"]
		tyrRInactiveIndex = np.array([bulkMoleculeIds.index(x) for x in tyrRInactiveId])
		tyrRInactiveCounts = bulkMoleculesReader.readColumn("counts")[:, tyrRInactiveIndex].reshape(-1)

		# Get the promoter-bound status of the tyrA gene
		tyrATfBoundId = ["EG11039_RNA__MONOMER0-162"]
		tyrATfBoundIndex = np.array([bulkMoleculeIds.index(x) for x in tyrATfBoundId])
		tyrATfBoundCounts = bulkMoleculesReader.readColumn("counts")[:, tyrATfBoundIndex].reshape(-1)

		# Get the amount of monomeric tyrA
		tyrAProteinId = ["CHORISMUTPREPHENDEHYDROG-MONOMER[c]"]
		tyrAProteinIndex = np.array([bulkMoleculeIds.index(x) for x in tyrAProteinId])
		tyrAProteinCounts = bulkMoleculesReader.readColumn("counts")[:, tyrAProteinIndex].reshape(-1)

		tyrAComplexId = ["CHORISMUTPREPHENDEHYDROG-CPLX[c]"]
		tyrAComplexIndex = np.array([bulkMoleculeIds.index(x) for x in tyrAComplexId])
		tyrAComplexCounts = bulkMoleculesReader.readColumn("counts")[:, tyrAComplexIndex].reshape(-1)

		bulkMoleculesReader.close()

		tyrAProteinTotalCounts = tyrAProteinCounts + 2 * tyrAComplexCounts

		# Compute the tyrA mass in the cell
		tyrAMw = sim_data.getter.getMass(tyrAProteinId)
		tyrAMass = 1. / nAvogadro * tyrAProteinTotalCounts * tyrAMw


		# Compute the proteome mass fraction
		proteomeMassFraction = tyrAMass.asNumber(units.fg) / proteinMass.asNumber(units.fg)

		# Get the tyrA synthesis probability
		rnaSynthProbReader = TableReader(os.path.join(simOutDir, "RnaSynthProb"))

		rnaIds = rnaSynthProbReader.readAttribute("rnaIds")
		tyrASynthProbId = ["EG11039_RNA[c]"]
		tyrASynthProbIndex = np.array([rnaIds.index(x) for x in tyrASynthProbId])
		tyrASynthProb = rnaSynthProbReader.readColumn("rnaSynthProb")[:, tyrASynthProbIndex].reshape(-1)

		tyrRBound = rnaSynthProbReader.readColumn("nActualBound")[:,tyrRIndex]
		
		rnaSynthProbReader.close()

		# Calculate total tyrR - active, inactive and bound
		tyrRTotalCounts = tyrRActiveCounts + tyrRInactiveCounts + tyrRBound

		# Compute moving averages
		width = 100

		tyrATfBoundCountsMA = np.convolve(tyrATfBoundCounts, np.ones(width) / width, mode = "same")
		tyrASynthProbMA = np.convolve(tyrASynthProb, np.ones(width) / width, mode = "same")
		
		##############################################################
		ax = plt.subplot(6, 1, 1)
		ax.plot(time, pheConcentration.asNumber(units.umol / units.L), color = "b")
		plt.ylabel("Internal PHE Conc. [uM]", fontsize = 6)

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
		ax.plot(time, tyrRActiveCounts)
		ax.plot(time, tyrRInactiveCounts)
		ax.plot(time, tyrRTotalCounts)
		plt.ylabel("TyrR Counts", fontsize = 6)
		plt.legend(["Active", "Inactive", "Total"], fontsize = 6)

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
		ax.plot(time, tyrATfBoundCounts, color = "b")
		ax.plot(time, tyrATfBoundCountsMA, color = "g")
		plt.ylabel("TyrR Bound To tyrA Promoter", fontsize = 6)

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
		ax.plot(time, tyrASynthProb, color = "b")
		ax.plot(time, tyrASynthProbMA, color = "g")
		plt.ylabel("tyrA Synthesis Prob.", fontsize = 6)

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
		ax.plot(time, tyrAProteinTotalCounts, color = "b")
		plt.ylabel("TyrA Counts", fontsize = 6)

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
		plt.ylabel("TyrA Mass Fraction of Proteome", fontsize = 6)

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
