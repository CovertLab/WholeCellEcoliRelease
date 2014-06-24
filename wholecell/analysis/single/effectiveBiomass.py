#!/usr/bin/env python
"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import division

import argparse
import os
import cPickle

import tables
import numpy as np
import matplotlib
# matplotlib.use("Agg")
from matplotlib import pyplot as plt

import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	with tables.open_file(os.path.join(simOutDir, "EffectiveBiomassObjective.hdf")) as h5file:
		timeStep = h5file.root.EffectiveBiomassObjective.col("time")
		effectiveBiomass = h5file.root.EffectiveBiomassObjective.col("effectiveBiomassObjective")
		nominalBiomass = h5file.root.EffectiveBiomassObjective.attrs.wildtypeBiomass
		metaboliteIds = h5file.root.EffectiveBiomassObjective.attrs.metaboliteIds

	ntpIds = ["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]
	ntpIdxs = np.array([np.where(ntpId == metaboliteIds)[0][0] for ntpId in ntpIds])

	aaIds = cPickle.load(open(kbFile)).aaIDs
	aaIdxs = np.array([np.where(aaId == metaboliteIds)[0][0] for aaId in aaIds])

	ntpEffective = effectiveBiomass[:, ntpIdxs].sum(1)

	aaEffective = effectiveBiomass[:, aaIdxs].sum(1)

	# normalizedEffective = np.nan_to_num(effectiveBiomass / effectiveBiomass.sum(1)[:, np.newaxis])
	# normalizedNominal = np.nan_to_num(nominalBiomass / nominalBiomass.sum())

	# relativeEffective =  np.nan_to_num(effectiveBiomass / nominalBiomass)
	# relNormEffective =  np.nan_to_num(normalizedEffective / normalizedNominal)

	# plt.figure(figsize = (8.5, 11))

	import ipdb; ipdb.set_trace()

	# plt.savefig(os.path.join(plotOutDir, plotOutFileName))


if __name__ == "__main__":
	# defaultKBFile = os.path.join(
	# 		wholecell.utils.constants.SERIALIZED_KB_DIR,
	# 		wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
	# 		)

	# parser = argparse.ArgumentParser()
	# parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	# parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	# parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	# parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	# args = parser.parse_args().__dict__

	# main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])

	# TODO: write this to provide meaningful output
	pass
