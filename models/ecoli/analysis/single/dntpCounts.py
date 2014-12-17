#!/usr/bin/env python
"""
Plot NTP counts

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf"))

	names = h.root.names

	moleculeIds = names.moleculeIDs.read()

	DNTP_IDS = ['DATP[c]', 'DCTP[c]', 'DGTP[c]', 'DTTP[c]']

	dntpIndexes = np.array([moleculeIds.index(dntpId) for dntpId in DNTP_IDS], np.int)
	bulkMolecules = h.root.BulkMolecules
	dntpCounts = bulkMolecules.read(0, None, 1, "counts")[:, dntpIndexes]


	DNMP_IDS = ['DAMP[n]', 'DCMP[n]', 'DGMP[n]', 'DTMP[n]']
	dnmpIndexes = np.array([moleculeIds.index(dntpId) for dntpId in DNMP_IDS], np.int)
	bulkMolecules = h.root.BulkMolecules
	dnmpCounts = bulkMolecules.read(0, None, 1, "counts")[:, dnmpIndexes]

	h.close()

	h = tables.open_file(os.path.join(simOutDir, "Mass.hdf"))
	table = h.root.Mass
	time = np.array([x["time"] for x in table.iterrows()])
	h.close()

	plt.figure(figsize = (8.5, 11))

	for idx in xrange(4):

		plt.subplot(2, 2, idx + 1)

		plt.plot(time / 60., dntpCounts[:, idx], linewidth = 2)
		plt.xlabel("Time (min)")
		plt.ylabel("Counts")
		plt.title(DNTP_IDS[idx])

		# print float(dntpCounts[-1, idx] + dnmpCounts[-1, idx]) / (dntpCounts[0, idx] + dnmpCounts[0, idx])
		# print float(dntpCounts[-1, idx]) / (dntpCounts[0, idx])


	plt.subplots_adjust(hspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)

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
