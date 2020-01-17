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

(sorted from sim_data.transcription.process.rnaData, mrna only, elements in range -1505:-1495)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/29/2015
"""

from __future__ import absolute_import

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
		all_mRNA_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_ids'))}

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

		rnaIndexes = np.array([all_mRNA_idx[x] for x in rnaIds], np.int)
		rnaCounts = mRNA_counts[:, rnaIndexes]

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		plt.figure(figsize = (8.5, 11))

		for subplotIdx in xrange(1, 10):

			plt.subplot(3, 3, subplotIdx)

			plt.plot(time / 60., rnaCounts[:, subplotIdx])
			plt.xlabel("Time (min)")
			plt.ylabel("mRNA counts")
			plt.title(names[subplotIdx].split(" - ")[0])

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
