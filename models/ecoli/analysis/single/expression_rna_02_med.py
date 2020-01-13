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

(sorted from sim_data.transcription.process.rnaData, mrna only, elements in range -125:-115)

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
		all_mRNA_ids = mRNA_counts_reader.readAttribute('mRNA_ids')

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

		rnaIndexes = np.array([all_mRNA_ids.index(x) for x in rnaIds], np.int)
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
