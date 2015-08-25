#!/usr/bin/env python
"""
Plots rnap capacity

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/18/2015
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
from wholecell.utils.sparkline import sparklineAxis, setAxisMaxMinY

FONT = {
		'size'	:	8
		}

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	nAvogadro = kb.constants.nAvogadro

	rnapSubunitIds = kb.process.complexation.getMonomers("APORNAP-CPLX[c]")['subunitIds']
	rnapSubunitStoich = kb.process.complexation.getMonomers("APORNAP-CPLX[c]")['subunitStoich']
	massFullRnapComplex = kb.getter.getMass(["APORNAP-CPLX[c]"])[0]
	rnapSubunitMasses = kb.getter.getMass(rnapSubunitIds)

	elongationRate = float(kb.growthRateParameters.rnaPolymeraseElongationRate.asNumber(units.nt / units.s))

	# Load rnap data
	rnapDataFile = TableReader(os.path.join(simOutDir, "RnapData"))

	actualElongations = rnapDataFile.readColumn("actualElongations")
	expectedElongations_recorded = rnapDataFile.readColumn("expectedElongations")
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	rnapDataFile.close()

	# Load count data for s30 proteins, rRNA, and final 30S complex
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	# Get indexes
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	rnapSubunitIndexes = np.array([moleculeIds.index(comp) for comp in rnapSubunitIds], np.int)

	# Load data
	rnapSubunitCounts = bulkMolecules.readColumn("counts")[:, rnapSubunitIndexes]

	bulkMolecules.close()

	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

	rnapIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRnaPoly")
	activeRnap = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, rnapIndex]

	uniqueMoleculeCounts.close()

	# Calculate statistics
	totalRnap = (activeRnap + np.min(rnapSubunitCounts / rnapSubunitStoich))
	totalRnapCapacity = totalRnap * elongationRate

	freeSubunitMass = (rnapSubunitMasses * rnapSubunitCounts / nAvogadro).asNumber(units.fg).sum(axis = 1)
	activeRnapMass = (massFullRnapComplex * activeRnap / nAvogadro).asNumber(units.fg)
	totalRnapMass = freeSubunitMass + activeRnapMass
	massFractionActive = activeRnapMass / totalRnapMass

	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	rnapCapacity_axis = plt.subplot(4,1,1)
	rnapCapacity_axis.plot(time / 60., totalRnapCapacity, label="Theoretical total rnap capacity", linewidth=2, color='b')
	rnapCapacity_axis.plot(time / 60., actualElongations, label="Actual elongations", linewidth=2, color='r')
	rnapCapacity_axis.set_ylabel("Amino acids polymerized")
	rnapCapacity_axis.legend(ncol=2)

	fractionalCapacity_axis = plt.subplot(4,1,2)
	fractionalCapacity_axis.plot(time / 60., actualElongations / totalRnapCapacity, label="Fraction of rnap capacity used", linewidth=2, color='k')
	fractionalCapacity_axis.set_ylabel("Fraction of rnap capacity used")
	fractionalCapacity_axis.set_yticks(np.arange(0., 1.05, 0.05))
	#fractionalCapacity_axis.get_yaxis().grid(b=True, which='major', color='b', linestyle='--')
	fractionalCapacity_axis.grid(b=True, which='major', color='b', linestyle='--')

	effectiveElongationRate_axis = plt.subplot(4,1,3)
	effectiveElongationRate_axis.plot(time / 60., actualElongations / activeRnap, label="Effective elongation rate", linewidth=2, color='k')
	effectiveElongationRate_axis.set_ylabel("Effective elongation rate (aa/s/rnap)")

	fractionActive_axis = plt.subplot(4,1,4)
	fractionActive_axis.plot(time / 60., massFractionActive, label="Mass fraction active", linewidth=2, color='k')
	fractionActive_axis.set_ylabel("Mass fraction of active rnaps")
	fractionActive_axis.set_yticks(np.arange(0., 1.1, 0.1))

	# Save
	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")

if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
