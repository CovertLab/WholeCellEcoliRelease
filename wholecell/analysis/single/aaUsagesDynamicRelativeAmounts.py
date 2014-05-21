#!/usr/bin/env python
"""
Plot NTP usages

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/2014
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

	h = tables.open_file(os.path.join(simOutDir, "AAUsage.hdf"))

	metaboliteIds = h.root.AAUsage._v_attrs["metaboliteIds"]
	normAAProductionBiomass = h.root.AAUsage._v_attrs["relativeAAProductionBiomass"]

	aaUsage = np.array(
		[x['translationAAUsageCurrent'] for x in h.root.AAUsage.iterrows()]
		)[1:, :]	# Ignore time point 0
	t = np.array(
		[x["time"] for x in h.root.AAUsage.iterrows()]
		)[1: ]	# Ignore time point 0

	h.close()

	normUsage = aaUsage / np.tile(
		aaUsage.sum(axis = 1).astype("float64").reshape(-1, 1), (1, 20)
		)

	plt.figure(figsize = (8.5, 11))

	for idx in xrange(20):

		plt.subplot(5, 4, idx + 1)

		plt.plot(t / 60., normUsage[:, idx], linewidth = 2)
		plt.plot([t[0] / 60., t[-1] / 60.], [normAAProductionBiomass[idx], normAAProductionBiomass[idx]], "r--")
		plt.xlabel("Time (min)")
		plt.ylabel("Relative usage")
		plt.title(metaboliteIds[idx])

	plt.subplots_adjust(hspace = 0.5)

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])