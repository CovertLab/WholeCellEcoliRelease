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

def main(simOutDir, plotOutDir, plotOutFileName):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf"))

	names = h.root.names

	moleculeIds = names.moleculeIDs.read()

	NTP_IDS = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']
	ntpIndexes = np.array([moleculeIds.index(ntpId) for ntpId in NTP_IDS], np.int)
	bulkMolecules = h.root.BulkMolecules
	ntpCounts = bulkMolecules.read(0, None, 1, "counts")[:, ntpIndexes]

	h.close()

	h = tables.open_file(os.path.join(simOutDir, "Mass.hdf"))
	table = h.root.Mass
	time = np.array([x["time"] for x in table.iterrows()])
	h.close()

	plt.figure(figsize = (8.5, 11))

	for idx in xrange(4):

		plt.subplot(2, 2, idx + 1)

		plt.plot(time / 60., ntpCounts[:, idx], linewidth = 2)
		plt.xlabel("Time (min)")
		plt.ylabel("Counts")
		plt.title(NTP_IDS[idx])

	plt.subplots_adjust(hspace = 0.5)

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])