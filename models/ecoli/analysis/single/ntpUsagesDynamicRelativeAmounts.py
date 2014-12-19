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

import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "NtpUsage.hdf"))

	metaboliteIds = h.root.NtpUsage._v_attrs["metaboliteIds"]
	normNtpProductionBiomass = h.root.NtpUsage._v_attrs["relativeNtpProductionBiomass"]
	relativeNtpUsage = h.root.NtpUsage._v_attrs["relativeNtpUsage"]

	ntpUsage = np.array(
		[x['transcriptionNtpUsageCurrent'] for x in h.root.NtpUsage.iterrows()]
		)[1:, :]	# Ignore time point 0
	t = np.array(
		[x["time"] for x in h.root.NtpUsage.iterrows()]
		)[1: ]	# Ignore time point 0

	h.close()

	normUsage = ntpUsage / np.tile(
		ntpUsage.sum(axis = 1).astype("float64").reshape(-1, 1), (1, 4)
		)

	plt.figure(figsize = (8.5, 11))

	for idx in xrange(4):

		plt.subplot(2, 2, idx + 1)

		plt.plot(t / 60., normUsage[:, idx], 'k', linewidth = 1)
		(biomass,) = plt.plot([t[0] / 60., t[-1] / 60.], [normNtpProductionBiomass[idx], normNtpProductionBiomass[idx]], "r--", linewidth = 2)
		(usages,) = plt.plot([t[0] / 60., t[-1] / 60.], [relativeNtpUsage[idx], relativeNtpUsage[idx]], "b--", linewidth = 2)
		plt.xlabel("Time (min)")
		plt.ylabel("Relative usage")
		plt.title(metaboliteIds[idx])

		if idx == 0:
			plt.legend([biomass, usages], ["Biomass production", "Expected usage"], loc = "best")

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
