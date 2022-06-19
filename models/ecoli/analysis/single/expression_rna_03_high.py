"""
Plot dynamic traces of genes with high expression (> 20 counts of mRNA)

EG10367_RNA	24.8	gapA	Glyceraldehyde 3-phosphate dehydrogenase
EG11036_RNA	25.2	tufA	Elongation factor Tu
EG50002_RNA	26.2	rpmA	50S Ribosomal subunit protein L27
EG10671_RNA	30.1	ompF	Outer membrane protein F
EG50003_RNA	38.7	acpP	Apo-[acyl carrier protein]
EG10669_RNA	41.1	ompA	Outer membrane protein A
EG10873_RNA	44.7	rplL	50S Ribosomal subunit protein L7/L12 dimer
EG12179_RNA	46.2	cspE	Transcription antiterminator and regulator of RNA stability
EG10321_RNA	53.2	fliC	Flagellin
EG10544_RNA	97.5	lpp		Murein lipoprotein
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
		all_mRNA_cistron_idx = {rna: i for i, rna in enumerate(mRNA_counts_reader.readAttribute('mRNA_cistron_ids'))}

		rnaIds = [
			"EG10367_RNA", "EG11036_RNA", "EG50002_RNA", "EG10671_RNA", "EG50003_RNA",
			"EG10669_RNA", "EG10873_RNA", "EG12179_RNA", "EG10321_RNA", "EG10544_RNA",
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

		cistron_indexes = np.array([all_mRNA_cistron_idx[x] for x in rnaIds], int)
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
