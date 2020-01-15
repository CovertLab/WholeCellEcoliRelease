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
