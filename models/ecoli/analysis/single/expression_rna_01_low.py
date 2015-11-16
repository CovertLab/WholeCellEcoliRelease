#!/usr/bin/env python
"""
Plot dynamic traces of genes with low expression (on average less than 1 mRNA)

G7355_RNA[c]	0.23	ypjD	Predicted inner membrane protein
EG11783_RNA[c]	0.23	intA	CP4-57 prophage; integrase
G7742_RNA[c]	0.23	yrfG	Purine nucleotidase
G6253_RNA[c]	0.23	ylaC	Predicted inner membrane protein
EG10632_RNA[c]	0.23	nagA	N-acetylglucosamine-6-phosphate deacetylase
EG11484_RNA[c]	0.23	yigZ	Predicted elongation factor
G7889_RNA[c]	0.23	lptG	LptG (part of LPS transport system)
EG10997_RNA[c]	0.24	mnmE	GTPase, involved in modification of U34 in tRNA
EG10780_RNA[c]	0.24	pspE	Thiosulfate sulfurtransferase
EG11060_RNA[c]	0.24	ushA	UDP-sugar hydrolase / 5'-ribonucleotidase / 5'-deoxyribonucleotidase

(sorted from kb.transcription.process.rnaData, mrna only, elements in range -1505:-1495)

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
		"G7355_RNA[c]", "EG11783_RNA[c]", "G7742_RNA[c]", "G6253_RNA[c]", "EG10632_RNA[c]",
		"EG11484_RNA[c]", "G7889_RNA[c]", "EG10997_RNA[c]", "EG10780_RNA[c]", "EG11060_RNA[c]",
		]
	names = [
		"ypjD - Predicted inner membrane protein",
		"intA - CP4-57 prophage; integrase",
		"yrfG - Purine nucleotidase",
		"ylaC - Predicted inner membrane protein",
		"nagA - N-acetylglucosamine-6-phosphate deacetylase",
		"yigZ - Predicted elongation factor",
		"lptG - LptG (part of LPS transport system)",
		"mnmE - GTPase, involved in modification of U34 in tRNA",
		"pspE - Thiosulfate sulfurtransferase",
		"ushA - UDP-sugar hydrolase / 5'-ribonucleotidase / 5'-deoxyribonucleotidase",
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
