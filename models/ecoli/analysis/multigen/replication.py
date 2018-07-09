from __future__ import absolute_import

import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

PLACE_HOLDER = -1

CRITICAL_N = [1, 2, 4, 8]


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(seedOutDir):
			raise Exception, "seedOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		sim_data = cPickle.load(open(simDataFile, "rb"))
		oriC = sim_data.constants.oriCCenter.asNumber()
		terC = sim_data.constants.terCCenter.asNumber()
		genomeLength = len(sim_data.process.replication.genome_sequence)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		# Get all cells
		allDir = ap.get_cells()

		fig, axesList = plt.subplots(6, sharex = True)
		fig.set_size_inches(11, 11)
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") / 60. / 60.
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			division_time = time.max()

			# Plot dna polymerase position
			sequenceIdx = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceIdx")
			sequenceLength = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("sequenceLength")
			reverseIdx = 1
			reverseCompIdx = 3
			reverseSequences = np.logical_or(sequenceIdx == reverseIdx, sequenceIdx == reverseCompIdx)
			sequenceLength[reverseSequences] = -1 * sequenceLength[reverseSequences]
			sequenceLength[sequenceLength == PLACE_HOLDER] = np.nan

			# Down sample dna polymerase position, every position is only plotted once here
			# using numpy ninja-ness
			unique, index, value = np.unique(sequenceLength, return_index=True, return_inverse=True)
			m = np.zeros_like(value, dtype=bool)
			m[index] = True
			m = m.reshape(sequenceLength.shape)
			sequenceLength[~m] = np.nan

			axesList[0].plot(time, sequenceLength, marker=',', markersize=1, linewidth=0)
			axesList[0].set_xlim([0, time.max()])
			axesList[0].set_yticks([-1 * genomeLength / 2, 0, genomeLength / 2])
			axesList[0].set_yticklabels(['-terC', 'oriC', '+terC'])
			axesList[0].set_ylabel("DNA polymerase\nposition (nt)")

			# Plot dna polymerase counts
			pairsOfForks = (sequenceIdx != PLACE_HOLDER).sum(axis = 1) / 4

			axesList[1].plot(time, pairsOfForks, linewidth=2)
			axesList[1].set_yticks(np.arange(0,7))
			axesList[1].set_ylim([0, 6])
			axesList[1].set_xlim([0, time.max()])
			axesList[1].plot([time.max(), time.max()], axesList[1].get_ylim(), 'k')
			axesList[1].set_ylabel("Pairs of\nforks")

			# Factors of critical initiation mass
			totalMass = mass.readColumn("cellMass")
			criticalInitiationMass = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalInitiationMass")
			criticalMassEquivalents = totalMass / criticalInitiationMass

			axesList[2].plot(time, criticalMassEquivalents, linewidth=2)
			for N in CRITICAL_N:
				axesList[2].plot([0, time.max()], [N]*2, 'k', linestyle='--')
			axesList[2].set_ylabel("Factors of critical\ninitiation mass")
			axesList[2].set_xlim([0, time.max()])

			# Dry mass
			dryMass = mass.readColumn("dryMass")
			axesList[3].plot(time, dryMass, linewidth = 2)
			axesList[3].set_ylabel("Dry mass (fg)")
			axesList[3].set_xlim([0, time.max()])
			mass.close()

			# Number of oriC
			numberOfOric = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("numberOfOric")
			axesList[4].plot(time, numberOfOric, linewidth=2)
			axesList[4].set_ylabel("Number of\noriC")
			axesList[4].set_xlim([0, time.max()])

			# Mass per oriC
			criticalMassPerOriC = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")
			axesList[5].plot(time, criticalMassPerOriC, linewidth=2)
			axesList[5].set_ylabel("Critical mass\nper oriC")
			axesList[5].set_xlim([0, time.max()])

		axesList[0].set_title("Replication plots")
		axesList[-1].set_xlabel("Time (hr)")

		plt.subplots_adjust(hspace = 0.2, wspace = 0.5)

		exportFigure(plt, plotOutDir, plotOutFileName,metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
