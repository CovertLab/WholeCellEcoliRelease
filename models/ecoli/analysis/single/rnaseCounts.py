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
import scipy.fftpack
import cPickle

import wholecell.utils.constants

from wholecell.io.tablereader import TableReader

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("moleculeIDs")

	RNase_IDS = ["EG10856-MONOMER[p]", "EG11620-MONOMER[c]", "EG10857-MONOMER[c]", "G7175-MONOMER[c]", "EG10858-MONOMER[c]", "EG10859-MONOMER[c]", "EG11299-MONOMER[c]", "EG10860-MONOMER[c]", "EG10861-MONOMER[c]", "G7365-MONOMER[c]", "EG10862-MONOMER[c]", "EG10863-MONOMER[c]", "EG11259-MONOMER[c]", "EG11547-MONOMER[c]", "EG10746-MONOMER[c]", "G7842-MONOMER[c]", "EG10743-MONOMER[c]", "EG10743_RNA[c]"]

	kb = cPickle.load(open(kbFile, "rb"))

	rnapRnaIndexes = np.array([moleculeIds.index(rnapRnaId) for rnapRnaId in RNase_IDS], np.int)
	rnapRnaCounts = bulkMolecules.readColumn("counts")[:, rnapRnaIndexes]
	time = bulkMolecules.readColumn("time")
	bulkMolecules.close()

	plt.figure(figsize = (8.5, 11))
	plt.rc('xtick', labelsize=5) 
	plt.rc('ytick', labelsize=5)

	for subplotIdx in xrange(0, 18):
		rnapRnaCountsIdx = subplotIdx
	
		plt.subplot(18, 1, 1 + subplotIdx)

		plt.plot(time / 60., rnapRnaCounts[:, rnapRnaCountsIdx])
		plt.xlabel("Time (min)", fontsize = 5)
		plt.ylabel("Protein counts", fontsize = 5)
		plt.title(RNase_IDS[rnapRnaCountsIdx], fontsize = 5)
		

		# identifying periodicity on RNA and protein copy numbers
		signal = rnapRnaCounts[:, rnapRnaCountsIdx]
		fourier = np.fft.fft(signal)

		# copmuting mean and std values of the power spectral density
		M = np.mean(abs(fourier))
		S = np.std(abs(fourier))

		# computing frequencies 
		n = signal.size
		timestep = 1 # second
		freq = np.fft.fftfreq(n, d=timestep)
		fft_freq = sorted(zip(abs(fourier),freq))[n-6:n]
		# import ipdb; ipdb.set_trace()

		# identifing peaks (frequency and period) with maximum PSD
		for i in xrange(0,len(fft_freq)):
			if fft_freq[i][1] > 0.: # only positive frequencies
				if 1. / fft_freq[i][1] < 3600.: # only periods lower than the doubling time
					if abs(fft_freq[i][0] - M) / S > 5: # strong and significant fft
						print RNase_IDS[rnapRnaCountsIdx], 1. / fft_freq[i][1] / 60.


	plt.subplots_adjust(hspace = 1.2, top = 0.95, bottom = 0.05)
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
