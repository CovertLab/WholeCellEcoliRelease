"""
Plots limited enzyme fluxes, protein counts, and transcription initiation events.

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/3/2017
"""

from __future__ import absolute_import

import os
import cPickle

import numpy as np
import matplotlib.pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from models.ecoli.processes.metabolism import COUNTS_UNITS, TIME_UNITS, VOLUME_UNITS
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()

		sim_data = cPickle.load(open(simDataFile, "rb"))
		enzymeComplexId = "CPLX0-8098[c]"
		enzymeMonomerId = "UGD-MONOMER[c]"
		enzymeRnaId = "G7091_RNA[c]"
		reactionId = "UGD-RXN"
		transcriptionFreq = 0.64
		metaboliteId = "UDP-GLUCURONATE[c]"

		simOutDir = os.path.join(allDir[0], "simOut")
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		enzymeComplexIndex = moleculeIds.index(enzymeComplexId)
		enzymeMonomerIndex = moleculeIds.index(enzymeMonomerId)
		metaboliteIndex = moleculeIds.index(metaboliteId)

		mRNA_counts_reader = TableReader(
			os.path.join(simOutDir, 'mRNACounts'))
		all_mRNA_ids = mRNA_counts_reader.readAttribute('mRNA_ids')
		enzymeRnaIndex = all_mRNA_ids.index(enzymeRnaId)

		rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
		rnap_data_rna_ids = rnapDataReader.readAttribute('rnaIds')
		enzyme_RNA_index_rnap_data = rnap_data_rna_ids.index(enzymeRnaId)

		time = []
		enzymeFluxes = []
		enzymeComplexCounts = []
		enzymeMonomerCounts = []
		enzymeRnaCounts = []
		enzymeRnaInitEvent = []
		metaboliteCounts = []

		for gen, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			time += TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeCounts = bulkMolecules.readColumn("counts")
			enzymeComplexCounts += moleculeCounts[:, enzymeComplexIndex].tolist()
			enzymeMonomerCounts += moleculeCounts[:, enzymeMonomerIndex].tolist()
			metaboliteCounts += moleculeCounts[:, metaboliteIndex].tolist()

			mRNA_counts_reader = TableReader(
				os.path.join(simOutDir, 'mRNACounts'))
			mRNA_counts = mRNA_counts_reader.readColumn('mRNA_counts')
			enzymeRnaCounts += mRNA_counts[:, enzymeRnaIndex].tolist()

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
			reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
			enzymeFluxes += reactionFluxes[:, np.where(reactionIDs == reactionId)[0][0]].tolist()

			rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
			enzymeRnaInitEvent += rnapDataReader.readColumn("rnaInitEvent")[:, enzyme_RNA_index_rnap_data].tolist()

		time = np.array(time)

		# Plot
		fig = plt.figure(figsize = (10, 10))
		rnaInitAxis = plt.subplot(6, 1, 1)
		rnaAxis = plt.subplot(6, 1, 2, sharex = rnaInitAxis)
		monomerAxis = plt.subplot(6, 1, 3, sharex = rnaInitAxis)
		complexAxis = plt.subplot(6, 1, 4, sharex = rnaInitAxis)
		fluxAxis = plt.subplot(6, 1, 5, sharex = rnaInitAxis)
		metAxis = plt.subplot(6, 1, 6, sharex = rnaInitAxis)

		rnaInitAxis.plot(time / 3600., enzymeRnaInitEvent)
		rnaInitAxis.set_title("%s transcription initiation events" % enzymeRnaId, fontsize = 10)
		rnaInitAxis.set_xlim([0, time[-1] / 3600.])
		rnaInitAxis.set_ylim([0, rnaInitAxis.get_ylim()[1] * 1.1])

		rnaAxis.plot(time / 3600., enzymeRnaCounts)
		rnaAxis.set_title("%s counts" % enzymeRnaId, fontsize = 10)

		monomerAxis.plot(time / 3600., enzymeMonomerCounts)
		monomerAxis.set_title("%s counts" % enzymeMonomerId, fontsize = 10)

		complexAxis.plot(time / 3600., enzymeComplexCounts)
		complexAxis.set_title("%s counts" % enzymeComplexId, fontsize = 10)

		fluxAxis.plot(time / 3600., enzymeFluxes)
		fluxAxis.set_title("%s flux (%s / %s / %s)" % (reactionId, COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS), fontsize = 10)

		metAxis.plot(time / 3600., metaboliteCounts)
		metAxis.set_title("%s counts" % metaboliteId, fontsize = 10)
		metAxis.set_xlabel("Time (hour)\n(%s frequency of at least 1 transcription per generation)" % transcriptionFreq, fontsize = 10)

		plt.subplots_adjust(wspace = 0.4, hspace = 0.4) #, right = 0.83, bottom = 0.05, left = 0.07, top = 0.95)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
