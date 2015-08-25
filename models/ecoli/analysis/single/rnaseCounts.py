#!/usr/bin/env python
"""
Plot RNAse counts

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/14/2015
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import scipy.fftpack
import cPickle

import wholecell.utils.constants

from wholecell.io.tablereader import TableReader

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")

	kb = cPickle.load(open(kbFile, "rb"))

	endoRnaseIds = kb.moleculeGroups.endoRnaseIds
	exoRnaseIds = kb.moleculeGroups.exoRnaseIds
	RNase_IDS = np.concatenate((endoRnaseIds, exoRnaseIds))

	# add mRNA PNPase
	PNP_RNA = ["EG10743_RNA[c]"]
	RNase_IDS = np.concatenate((RNase_IDS, PNP_RNA))	

	rnapRnaIndexes = np.array([moleculeIds.index(rnapRnaId) for rnapRnaId in RNase_IDS], np.int)
	rnapRnaCounts = bulkMolecules.readColumn("counts")[:, rnapRnaIndexes]
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	bulkMolecules.close()

	plt.figure(figsize = (8.5, 11))
	plt.rc('xtick', labelsize=5) 
	plt.rc('ytick', labelsize=5)

	for subplotIdx in xrange(0, 19):
		rnapRnaCountsIdx = subplotIdx
	
		plt.subplot(19, 1, 1 + subplotIdx)

		plt.plot(time / 60., rnapRnaCounts[:, rnapRnaCountsIdx])
		plt.xlabel("Time (min)", fontsize = 5)
		plt.ylabel("Protein counts", fontsize = 5)
		plt.title(RNase_IDS[rnapRnaCountsIdx], fontsize = 5)

		signal = rnapRnaCounts[:, rnapRnaCountsIdx]
		if subplotIdx == 17:
			np.savetxt(os.path.join(plotOutDir, 'PNPase-MONOMER[c].txt'), signal)
		if subplotIdx == 18:
			np.savetxt(os.path.join(plotOutDir, 'PNPase-RNA[c].txt'), signal)

		# identifying periodicity on RNA and protein copy numbers
		fourier = np.fft.fft(signal)

		# copmuting mean and std values of the power spectral density
		M = np.mean(abs(fourier))
		S = np.std(abs(fourier))

		# computing frequencies 
		n = signal.size
		timestep = 1 # second
		freq = np.fft.fftfreq(n, d=timestep)
		fft_freq = sorted(zip(abs(fourier),freq))[n-6:n]

		# identifing peaks (frequency and period) with maximum PSD
		for i in xrange(0,len(fft_freq)):
			if fft_freq[i][1] > 0.: # only positive frequencies
				if 1. / fft_freq[i][1] < 3600.: # only periods lower than the doubling time
					if abs(fft_freq[i][0] - M) / S > 3: # strong and significant fft
						print RNase_IDS[rnapRnaCountsIdx], 1. / fft_freq[i][1] / 60. # period (min)


	plt.subplots_adjust(hspace = 1.2, top = 0.95, bottom = 0.05)
	plt.savefig(os.path.join(plotOutDir, plotOutFileName))


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
