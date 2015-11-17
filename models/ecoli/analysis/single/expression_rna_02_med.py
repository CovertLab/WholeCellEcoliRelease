#!/usr/bin/env python
"""
Plot dynamic traces of genes with medium expression (about 5 counts of mRNA)

EG10789_RNA[c]	4.9	ptsI	PTS enzyme I
EG11556_RNA[c]	5.0	talB	Transaldolase
EG12095_RNA[c]	5.0	secG	SecG
G1_RNA[c]		5.0	thiS	ThiS protein
G360_RNA[c]		5.1	flgD	Flagellar biosynthesis
EG10944_RNA[c]	5.1	serA	(S)-2-hydroxyglutarate reductase
EG12419_RNA[c]	5.2	gatY	GatY
EG10372_RNA[c]	5.2	gdhA	Glutamate dehydrogenase
EG10104_RNA[c]	5.2	atpG	ATP synthase F1 complex - gamma subunit
EG10539_RNA[c]	5.3	livJ	Branched chain amino acid ABC transporter - periplasmic binding protein

(sorted from kb.transcription.process.rnaData, mrna only, elements in range -125:-115)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/29/2015
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	rnaIds = [
		"EG10789_RNA[c]", "EG11556_RNA[c]", "EG12095_RNA[c]", "G1_RNA[c]", "G360_RNA[c]",
		"EG10944_RNA[c]", "EG12419_RNA[c]", "EG10372_RNA[c]", "EG10104_RNA[c]", "EG10539_RNA[c]",
		]
	names = [
		"ptsI - PTS enzyme I",
		"talB - Transaldolase",
		"secG - SecG",
		"thiS - ThiS protein",
		"flgD - Flagellar biosynthesis",
		"serA - (S)-2-hydroxyglutarate reductase",
		"gatY - GatY",
		"gdhA - Glutamate dehydrogenase",
		"atpG - ATP synthase F1 complex - gamma subunit",
		"livJ - Branched chain amino acid ABC transporter - periplasmic binding protein",
	]

	moleculeIds = bulkMolecules.readAttribute("objectNames")
	rnaIndexes = np.array([moleculeIds.index(x) for x in rnaIds], np.int)
	rnaCounts = bulkMolecules.readColumn("counts")[:, rnaIndexes]

	bulkMolecules.close()

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	plt.figure(figsize = (8.5, 11))

	for subplotIdx in xrange(1, 10):

		plt.subplot(3, 3, subplotIdx)

		plt.plot(time / 60., rnaCounts[:, subplotIdx])
		plt.xlabel("Time (min)")
		plt.ylabel("mRNA counts")
		plt.title(names[subplotIdx].split(" - ")[0])

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
