#!/usr/bin/env python
"""
Plot RNAse counts

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/14/2015
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import wholecell.utils.constants

from wholecell.containers.unique_molecules_data import UniqueMoleculesData

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	h = tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf"))

	names = h.root.names
	bulkMolecules = h.root.BulkMolecules

	moleculeIds = names.moleculeIDs.read()
	#rnapId = "EG12146-MONOMER[c]"
	#rnapIndex = moleculeIds.index(rnapId)
	#rnapCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, rnapIndex]

	RNAP_RNA_IDS = ["EG10856-MONOMER[p]", "EG11620-MONOMER[c]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10858-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]", "EG10863-MONOMER[c]", "EG11259-MONOMER[c]", "EG11547-MONOMER[c]"]
	rnapRnaIndexes = np.array([moleculeIds.index(rnapRnaId) for rnapRnaId in RNAP_RNA_IDS], np.int)
	rnapRnaCounts = bulkMolecules.read(0, None, 1, "counts")[:, rnapRnaIndexes]
	#import ipdb
	#ipdb.set_trace()

	h.close()

	h = tables.open_file(os.path.join(simOutDir, "UniqueMoleculeCounts.hdf"))

	uniqueMoleculeCounts = h.root.UniqueMoleculeCounts
	time = uniqueMoleculeCounts.col("time")

	h.close()

	plt.figure(figsize = (8.5, 11))
	plt.rc('xtick', labelsize=5) 
	plt.rc('ytick', labelsize=5)

	for subplotIdx in xrange(0, 14):
		rnapRnaCountsIdx = subplotIdx
	
		plt.subplot(14, 1, 1 + subplotIdx)

		plt.plot(time / 60., rnapRnaCounts[:, rnapRnaCountsIdx])
		plt.xlabel("Time (min)", fontsize = 5)
		plt.ylabel("Protein counts", fontsize = 5)
		plt.title(RNAP_RNA_IDS[rnapRnaCountsIdx], fontsize = 5)

	plt.subplots_adjust(hspace = 1, top = 0.95, bottom = 0.05)
	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

	# h.close()

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
