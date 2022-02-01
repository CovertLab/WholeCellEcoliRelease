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
from six.moves import cPickle, range

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot
from six.moves import zip


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile, "rb"))
		all_cistron_ids = sim_data.process.transcription.cistron_data["id"].tolist()

		cistron_ids = [
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

		cistron_idxs = [all_cistron_ids.index(x) for x in cistron_ids]
		deg_rates = sim_data.process.transcription.cistron_data['deg_rate'][cistron_idxs]


		# Get all cells
		allDir = self.ap.get_cells()

		rna_cistron_degraded_counts = []
		rna_cistron_counts = []
		dts = []

		N = 100
		for idx, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			main_reader = TableReader(os.path.join(simOutDir, "Main"))
			initialTime = main_reader.readAttribute("initialTime")
			time = main_reader.readColumn("time") - initialTime
			N = np.fmin(N, time.size)

			dts.append(main_reader.readColumn("timeStepSec"))

			rnaDegradationListener = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
			rna_cistron_degraded_counts.append(rnaDegradationListener.readColumn('count_RNA_degraded_per_cistron')[:, cistron_idxs])

			mRNA_counts_reader = TableReader(
				os.path.join(simOutDir, 'mRNACounts'))
			all_mRNA_cistron_ids = mRNA_counts_reader.readAttribute('mRNA_cistron_ids')
			cistron_indexes = np.array([all_mRNA_cistron_ids.index(x) for x in cistron_ids], int)
			rna_cistron_counts.append(mRNA_counts_reader.readColumn("mRNA_cistron_counts")[:, cistron_indexes])

		rnaDegradedCountsAveraged = []
		rnaCountsAveraged = []

		for dt, rnaDegradedCount, rnaCount in zip(dts, rna_cistron_degraded_counts, rna_cistron_counts):
			tmpArray = np.nan * np.ones_like(rnaDegradedCount)
			for colIdx in range(tmpArray.shape[1]):
				tmpArray[:, colIdx] = np.convolve(rnaDegradedCount[:, colIdx] / dt, np.ones(N) / N, mode = "same")
			rnaDegradedCountsAveraged.append(tmpArray[N:-1*N, :])

			tmpArray = np.nan * np.ones_like(rnaCount)
			for colIdx in range(tmpArray.shape[1]):
				tmpArray[:, colIdx] = np.convolve(rnaCount[:, colIdx], np.ones(N) / N, mode = "same")
			rnaCountsAveraged.append(tmpArray[N:-1*N, :])

		rnaDegradedCountsAveraged = np.vstack(rnaDegradedCountsAveraged)
		rnaCountsAveraged = np.vstack(rnaCountsAveraged)

		plt.figure(figsize = (8.5, 11))

		for subplotIdx in range(1, 10):

			plt.subplot(3, 3, subplotIdx)

			y = rnaDegradedCountsAveraged[:, subplotIdx]

			A = rnaCountsAveraged[:, subplotIdx]
			try:
				kdeg, _, _, _ = np.linalg.lstsq(A[:, np.newaxis], y, rcond=None)
			except ValueError:
				# TODO: Come up with a better/more descriptive error message
				# This is to handle errors that occurs when running short simulations
				print("Skipping subplot %d because not enough data" % (subplotIdx,))
				continue

			plt.scatter(
				rnaCountsAveraged[::N, subplotIdx],
				rnaDegradedCountsAveraged[::N, subplotIdx]
					)
			# plt.plot(time / 60, y)
			plt.xlabel("RNA (counts)", size = 10)
			plt.ylabel("RNA degraded (counts)", size = 10)
			plt.title(names[subplotIdx].split(" - ")[0] +
				"\n" +
				"kdeg meas: %0.1e\n" % (kdeg,) +
				"kdeg exp:  %0.1e" % deg_rates[subplotIdx].asNumber(1 / units.s),
				size = 10,
				)

		plt.subplots_adjust(hspace = 0.5, top = 0.95, bottom = 0.05)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
