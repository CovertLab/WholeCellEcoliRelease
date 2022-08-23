"""
Plot dynamic traces of genes with medium expression (about 5 counts of mRNA)

EG10789_RNA	4.9	ptsI	PTS enzyme I
EG11556_RNA	5.0	talB	Transaldolase
EG12095_RNA	5.0	secG	SecG
G1_RNA		5.0	thiS	ThiS protein
G360_RNA	5.1	flgD	Flagellar biosynthesis
EG10944_RNA	5.1	serA	(S)-2-hydroxyglutarate reductase
EG12419_RNA	5.2	gatY	GatY
EG10372_RNA	5.2	gdhA	Glutamate dehydrogenase
EG10104_RNA	5.2	atpG	ATP synthase F1 complex - gamma subunit
EG10539_RNA	5.3	livJ	Branched chain amino acid ABC transporter - periplasmic binding protein

(sorted from sim_data.transcription.process.rnaData, mrna only, elements in range -125:-115)
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from six.moves import range


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		mRNA_counts_reader = TableReader(os.path.join(simOutDir, 'mRNACounts'))
		mRNA_cistron_counts = mRNA_counts_reader.readColumn('mRNA_cistron_counts')
		all_mRNA_cistron_idx = {cistron: i for i, cistron in enumerate(mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}

		cistron_ids = [
			"EG10789_RNA", "EG11556_RNA", "EG12095_RNA", "G1_RNA", "G360_RNA",
			"EG10944_RNA", "EG12419_RNA", "EG10372_RNA", "EG10104_RNA", "EG10539_RNA",
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

		cistron_indexes = np.array([all_mRNA_cistron_idx[x] for x in cistron_ids], int)
		cistron_counts = mRNA_cistron_counts[:, cistron_indexes]

		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		plt.figure(figsize = (8.5, 11))

		for subplotIdx in range(1, 10):

			plt.subplot(3, 3, subplotIdx)

			plt.plot(time / 60., cistron_counts[:, subplotIdx])
			plt.xlabel("Time (min)")
			plt.ylabel("mRNA counts")
			plt.title(names[subplotIdx].split(" - ")[0])

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
