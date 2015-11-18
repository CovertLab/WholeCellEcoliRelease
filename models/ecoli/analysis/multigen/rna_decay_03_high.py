#!/usr/bin/env python
"""
Plot dynamic traces of genes with high expression (> 20 counts of mRNA)

EG10367_RNA[c]	24.8	gapA	Glyceraldehyde 3-phosphate dehydrogenase
EG11036_RNA[c]	25.2	tufA	Elongation factor Tu
EG50002_RNA[c]	26.2	rpmA	50S Ribosomal subunit protein L27
EG10671_RNA[c]	30.1	ompF	Outer membrane protein F
EG50003_RNA[c]	38.7	acpP	Apo-[acyl carrier protein]
EG10669_RNA[c]	41.1	ompA	Outer membrane protein A
EG10873_RNA[c]	44.7	rplL	50S Ribosomal subunit protein L7/L12 dimer
EG12179_RNA[c]	46.2	cspE	Transcription antiterminator and regulator of RNA stability
EG10321_RNA[c]	53.2	fliC	Flagellin
EG10544_RNA[c]	97.5	lpp		Murein lipoprotein

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/29/2015
"""

import argparse
import os

import numpy as np
import cPickle
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
import wholecell.utils.constants
from wholecell.utils import units

def main(seedOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(seedOutDir):
		raise Exception, "seedOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	kb = cPickle.load(open(kbFile, "rb"))
	allRnaIds = kb.process.transcription.rnaData["id"].tolist()

	rnaIds = [
		"EG10367_RNA[c]", "EG11036_RNA[c]", "EG50002_RNA[c]", "EG10671_RNA[c]", "EG50003_RNA[c]",
		"EG10669_RNA[c]", "EG10873_RNA[c]", "EG12179_RNA[c]", "EG10321_RNA[c]", "EG10544_RNA[c]",
		]
	names = [
		"gapA - Glyceraldehyde 3-phosphate dehydrogenase",
		"tufA - Elongation factor Tu",
		"rpmA - 50S Ribosomal subunit protein L27",
		"ompF - Outer membrane protein F",
		"acpP - Apo-[acyl carrier protein]",
		"ompA - Outer membrane protein A",
		"rplL - 50S Ribosomal subunit protein L7/L12 dimer",
		"cspE - Transcription antiterminator and regulator of RNA stability",
		"fliC - Flagellin",
		"lpp - Murein lipoprotein",
	]

	rnaIdxs = [allRnaIds.index(x) for x in rnaIds]
	degRates = kb.process.transcription.rnaData["degRate"][rnaIdxs]

	ap = AnalysisPaths(seedOutDir)

	# Get all cells
	allDir = ap.getAll()

	rnaDegradedCounts = []
	rnaCounts = []

	N = 100
	for idx, simDir in enumerate(allDir):
		simOutDir = os.path.join(simDir, "simOut")
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
		N = np.fmin(N, time.size)

		rnaDegradationListener = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
		rnaDegradedCounts.append(rnaDegradationListener.readColumn('countRnaDegraded')[:, rnaIdxs])
		rnaDegradationListener.close()

		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		rnaIndexes = np.array([moleculeIds.index(x) for x in rnaIds], np.int)
		rnaCounts.append(bulkMolecules.readColumn("counts")[:, rnaIndexes])
		bulkMolecules.close()

	rnaDegradedCountsAveraged = []
	rnaCountsAveraged = []

	for rnaDegradedCount, rnaCount in zip(rnaDegradedCounts, rnaCounts):
		tmpArray = np.nan * np.ones_like(rnaDegradedCount)
		for colIdx in xrange(tmpArray.shape[1]):
			tmpArray[:, colIdx] = np.convolve(rnaDegradedCount[:, colIdx], np.ones(N) / N, mode = "same")
		rnaDegradedCountsAveraged.append(tmpArray[N:-1*N, :])

		tmpArray = np.nan * np.ones_like(rnaCount)
		for colIdx in xrange(tmpArray.shape[1]):
			tmpArray[:, colIdx] = np.convolve(rnaCount[:, colIdx], np.ones(N) / N, mode = "same")
		rnaCountsAveraged.append(tmpArray[N:-1*N, :])

	rnaDegradedCountsAveraged = np.vstack(rnaDegradedCountsAveraged)
	rnaCountsAveraged = np.vstack(rnaCountsAveraged)

	plt.figure(figsize = (8.5, 11))

	for subplotIdx in xrange(1, 10):

		plt.subplot(3, 3, subplotIdx)

		y = rnaDegradedCountsAveraged[:, subplotIdx]

		A = rnaCountsAveraged[:, subplotIdx]
		kdeg, _, _, _ = np.linalg.lstsq(A[:, np.newaxis], y)

		plt.scatter(
			rnaCountsAveraged[::N, subplotIdx],
			rnaDegradedCountsAveraged[::N, subplotIdx]
				)
		# plt.plot(time / 60, y)
		plt.xlabel("RNA (counts)", size = 10)
		plt.ylabel("RNA degraded (counts)", size = 10)
		plt.title(names[subplotIdx].split(" - ")[0] +
			"\n" +
			"kdeg meas: %0.1e\n" % kdeg +
			"kdeg exp:  %0.1e" % degRates[subplotIdx].asNumber(1 / units.s),
			size = 10,
			)

	plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
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
