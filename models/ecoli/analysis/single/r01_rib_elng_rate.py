#!/usr/bin/env python
"""
Plots ribosome capacity

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/20/2014
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
	ribosomeSubunitIds = []
	ribosomeSubunitIds.extend(kb.moleculeGroups.s50_fullComplex)
	ribosomeSubunitIds.extend(kb.moleculeGroups.s30_fullComplex)
	ribosomeSubunitIds.extend(kb.moleculeGroups.s50_proteinComplexes)
	ribosomeSubunitIds.extend(kb.process.complexation.getMonomers(kb.moleculeGroups.s50_fullComplex[0])['subunitIds'])
	ribosomeSubunitIds.extend(kb.process.complexation.getMonomers(kb.moleculeGroups.s30_fullComplex[0])['subunitIds'])
	ribosomeSubunitMasses = kb.getter.getMass(ribosomeSubunitIds)
	mass70s = (kb.getter.getMass(kb.moleculeGroups.s50_fullComplex) + kb.getter.getMass(kb.moleculeGroups.s30_fullComplex))[0]

	elongationRate = float(kb.growthRateParameters.ribosomeElongationRate.asNumber(units.aa / units.s))

	# Load ribosome data
	massFile = TableReader(os.path.join(simOutDir, "RibosomeData"))

	actualElongations = massFile.readColumn("actualElongations")
	expectedElongations_recorded = massFile.readColumn("expectedElongations")
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	massFile.close()

	# Load count data for s30 proteins, rRNA, and final 30S complex
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	# Get indexes
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	ribosomeSubunitIndexes = np.array([moleculeIds.index(comp) for comp in ribosomeSubunitIds], np.int)

	# Load data
	ribosomeSubunitCounts = bulkMolecules.readColumn("counts")[:, ribosomeSubunitIndexes]

	bulkMolecules.close()

	uniqueMoleculeCounts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))

	ribosomeIndex = uniqueMoleculeCounts.readAttribute("uniqueMoleculeIds").index("activeRibosome")
	activeRibosome = uniqueMoleculeCounts.readColumn("uniqueMoleculeCounts")[:, ribosomeIndex]

	uniqueMoleculeCounts.close()

	# Calculate statistics
	totalRibosome = (activeRibosome + np.min(ribosomeSubunitCounts))
	totalRibosomeCapacity = totalRibosome * elongationRate

	freeSubunitMass = (ribosomeSubunitMasses * ribosomeSubunitCounts / nAvogadro).asNumber(units.fg).sum(axis = 1)
	activeRibosomeMass = (mass70s * activeRibosome / nAvogadro).asNumber(units.fg)
	totalRibosomeMass = freeSubunitMass + activeRibosomeMass
	massFractionActive = activeRibosomeMass / totalRibosomeMass

	f = plt.figure(figsize = (1.25, 1.))
	#f = plt.figure(figsize = (6.25, 6.))
	ax = f.add_axes([0, 0, 1, 1])
	ax.axis("off")

	eff_elng_rate = actualElongations / activeRibosome
	ax.plot(time, eff_elng_rate, label="Effective elongation rate", linewidth=1, color='k')
	ax.set_ylim([np.nanmin(eff_elng_rate) - 0.5, np.nanmax(eff_elng_rate) + 0.5])
	print np.nanmean(eff_elng_rate)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, "R01_effective_ribosome_elongation_rate")
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
