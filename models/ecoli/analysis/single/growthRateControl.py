#!/usr/bin/env python
"""
Plots ribosome capacity

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 12/14/2015
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

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	sim_data = cPickle.load(open(simDataFile, "rb"))
	max_elongationRate = sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s)
	elongationRate = float(sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))

	# Load ribosome data
	ribosomeDataFile = TableReader(os.path.join(simOutDir, "RibosomeData"))

	effectiveElongationRate = ribosomeDataFile.readColumn("effectiveElongationRate")
	rrnInitRate = ribosomeDataFile.readColumn("rrnInitRate")
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	ribosomeDataFile.close()

	# Load bulk molecule data
	bulkMoleculesFile = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	bulkMoleculeIds = bulkMoleculesFile.readAttribute("objectNames")
	rrn_idx = bulkMoleculesFile.readAttribute("objectNames").index('rrn_operon')
	rrn_counts = bulkMoleculesFile.readColumn("counts")[:, rrn_idx]
	bulkMoleculesFile.close()

	plt.figure(figsize = (8.5, 11))

	effectiveElongationRate_axis = plt.subplot(3,1,1)
	effectiveElongationRate_axis.plot(time / 60., effectiveElongationRate, label="Effective elongation rate", linewidth=2, color='k')
	effectiveElongationRate_axis.plot(time / 60., max_elongationRate * np.ones(time.size), 'r--')
	effectiveElongationRate_axis.set_ylabel("Effective elongation rate (aa/s/ribosome)")

	rrnCounts_axis = plt.subplot(3,1,2)
	rrnCounts_axis.plot(time / 60., rrn_counts, label="Rrn operon counts", linewidth=2, color='k')
	rrnCounts_axis.set_ylim([rrn_counts.min() - 1, rrn_counts.max() + 1])
	rrnCounts_axis.set_ylabel("Rrn operons")

	rrnCounts_axis = plt.subplot(3,1,3)
	rrnCounts_axis.plot(time / 60., rrnInitRate, label="Rrn init", linewidth=2, color='k')
	rrnCounts_axis.set_ylabel("Rrn init rate")

	# Save
	plt.subplots_adjust(hspace = 0.5, wspace = 0.6)

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
