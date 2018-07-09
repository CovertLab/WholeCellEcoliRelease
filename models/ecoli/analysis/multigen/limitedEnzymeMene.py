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

PLOT_DOWNSTREAM_ENZYME = False


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		return

		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		enzymeComplexId = "MENE-CPLX[c]"
		enzymeMonomerId = "O-SUCCINYLBENZOATE-COA-LIG-MONOMER[c]"
		enzymeRnaId = "EG11532_RNA[c]"
		reactionId = "O-SUCCINYLBENZOATE-COA-LIG-RXN"
		transcriptionFreq = 0.84
		metaboliteIds = ["REDUCED-MENAQUINONE[c]", "CPD-12115[c]"]
		downstreamEnzyme = "DMK-MONOMER[i]"

		# Get all cells
		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)
		allDir = ap.get_cells()
		if allDir.tolist() == []:
			return

		sim_data = cPickle.load(open(simDataFile, "rb"))
		rnaIds = sim_data.process.transcription.rnaData["id"]
		isMRna = sim_data.process.transcription.rnaData["isMRna"]
		mRnaIndexes = np.where(isMRna)[0]
		mRnaIds = np.array([rnaIds[x] for x in mRnaIndexes])

		simOutDir = os.path.join(allDir[0], "simOut")
		bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		moleculeIds = bulkMolecules.readAttribute("objectNames")
		enzymeComplexIndex = moleculeIds.index(enzymeComplexId)
		enzymeMonomerIndex = moleculeIds.index(enzymeMonomerId)
		enzymeRnaIndex = moleculeIds.index(enzymeRnaId)
		metaboliteIndexes = [moleculeIds.index(x) for x in metaboliteIds]
		downstreamEnzymeIndex = moleculeIds.index(downstreamEnzyme)
		bulkMolecules.close()

		time = []
		enzymeFluxes = []
		enzymeComplexCounts = []
		enzymeMonomerCounts = []
		enzymeRnaCounts = []
		enzymeRnaInitEvent = []
		metaboliteCounts = np.array([])
		downstream = []

		for gen, simDir in enumerate(allDir):
			simOutDir = os.path.join(simDir, "simOut")

			time += TableReader(os.path.join(simOutDir, "Main")).readColumn("time").tolist()

			bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
			moleculeCounts = bulkMolecules.readColumn("counts")
			enzymeComplexCounts += moleculeCounts[:, enzymeComplexIndex].tolist()
			enzymeMonomerCounts += moleculeCounts[:, enzymeMonomerIndex].tolist()
			enzymeRnaCounts += moleculeCounts[:, enzymeRnaIndex].tolist()
			if gen == 0:
				metaboliteCounts = moleculeCounts[:, metaboliteIndexes]
			else:
				metaboliteCounts = np.vstack((metaboliteCounts, moleculeCounts[:, metaboliteIndexes]))

			if PLOT_DOWNSTREAM_ENZYME:
				downstream += moleculeCounts[:, downstreamEnzymeIndex].tolist()
			bulkMolecules.close()

			fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
			reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
			reactionFluxes = np.array(fbaResults.readColumn("reactionFluxes"))
			enzymeFluxes += reactionFluxes[:, np.where(reactionIDs == reactionId)[0][0]].tolist()
			fbaResults.close()

			rnapDataReader = TableReader(os.path.join(simOutDir, "RnapData"))
			enzymeRnaInitEvent += rnapDataReader.readColumn("rnaInitEvent")[:, np.where(mRnaIds == enzymeRnaId)[0][0]].tolist()
			rnapDataReader.close()

		time = np.array(time)


		# Plot
		nrows = 6
		if PLOT_DOWNSTREAM_ENZYME:
			nrows = 7
		fig = plt.figure(figsize = (10, 10))
		rnaInitAxis = plt.subplot(nrows, 1, 1)
		rnaAxis = plt.subplot(nrows, 1, 2, sharex = rnaInitAxis)
		monomerAxis = plt.subplot(nrows, 1, 3, sharex = rnaInitAxis)
		complexAxis = plt.subplot(nrows, 1, 4, sharex = rnaInitAxis)
		fluxAxis = plt.subplot(nrows, 1, 5, sharex = rnaInitAxis)
		metAxis = plt.subplot(nrows, 1, 6, sharex = rnaInitAxis)

		if PLOT_DOWNSTREAM_ENZYME:
			axis = plt.subplot(nrows, 1, 7, sharex = rnaInitAxis)
			axis.plot(time / 3600., downstream)
			axis.set_title("%s counts" % downstreamEnzyme)

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
		metAxis.set_title("Metabolite counts (%s)" % metaboliteIds, fontsize = 10)
		metAxis.set_xlabel("Time (hour)\n(%s frequency of at least 1 transcription per generation)" % transcriptionFreq, fontsize = 10)

		plt.subplots_adjust(wspace = 0.4, hspace = 0.4)
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
