#!/usr/bin/env python
"""
Plots various effects that may be limiting growth

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

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	nAvogadro = kb.constants.nAvogadro

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Load data
	growthLimitsDataFile = TableReader(os.path.join(simOutDir, "GrowthLimits"))
	fractionAAsUsed = growthLimitsDataFile.readColumn("fractionAAsUsed")
	fractionGtpLimit = growthLimitsDataFile.readColumn("fractionGtpLimit")
	gtpPoolSize = growthLimitsDataFile.readColumn("gtpPoolSize")
	gtpRequestSize = growthLimitsDataFile.readColumn("gtpRequestSize")
	gtpAllocated = growthLimitsDataFile.readColumn("gtpAllocated")
	gtpPerElongation = growthLimitsDataFile.readColumn("gtpPerElongation")
	growthLimitsDataFile.close()


	# Calculate concentration data
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	bulkMoleculeIds = bulkMolecules.readAttribute("objectNames")

	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = units.fg * mass.readColumn("cellMass")

	cellDensity = kb.constants.cellDensity

	concentrationSetpoints = kb.process.metabolism.metabolitePoolConcentrations

	poolIds = np.array(kb.process.metabolism.metabolitePoolIDs)
	poolIndexes = np.array([bulkMoleculeIds.index(x) for x in poolIds])
	poolCounts = bulkMolecules.readColumn("counts")[:, poolIndexes]
	poolMols = 1/nAvogadro * poolCounts
	volume = cellMass / cellDensity
	poolConcentrations = poolMols * 1/volume[:,np.newaxis]

	# Compare
	common_units = units.mmol / units.L
	concSetpoint = np.tile(concentrationSetpoints.asNumber(common_units),(time.size,1))
	poolConc = poolConcentrations.asNumber(common_units)
	comp = np.isclose(concSetpoint, poolConc)
	concentrations_that_shift = np.logical_not(np.all(comp, axis = 0))

	bulkMolecules.close()


	fig = plt.figure(figsize = (8.5, 11))
	ax = plt.subplot(6,1,1)
	ax.plot(time / 60., fractionGtpLimit, linewidth=2, color='b')
	ax.set_ylabel("gtp used /\n gtp allocated")
	ax.set_title("Fraction of GTP limit used by translation")
	ax.text(0.5,0.5, "gtp per elongation {}".format(gtpPerElongation.max()))

	ax = plt.subplot(6,1,2)
	ax.plot(time / 60., fractionAAsUsed, linewidth=2, color='b')
	ax.set_ylabel("aa used /\n aa allocated")
	ax.set_title("Fraction of AA allocated used by translation")

	ax = plt.subplot(6,1,3)
	ax.plot(time / 60., gtpPoolSize, linewidth=2, label="pool size", color='k')
	ax.plot(time / 60., gtpRequestSize, linewidth=2, label="request size", color='b')
	ax.plot(time / 60., gtpAllocated, linewidth=2, label="allocated", color='r', linestyle='--')
	ax.legend()
	ax.set_ylabel("GTP")

	ax = plt.subplot(6,1,4)
	idx = 3
	ax.plot(time / 60., poolConc[:,idx], linewidth=2, label="pool size", color='k')
	ax.plot(time / 60., concSetpoint[:,idx], linewidth=2)
	ax.set_ylabel('{} mmol / L'.format(poolIds[idx]))

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
