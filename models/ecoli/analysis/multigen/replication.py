from __future__ import absolute_import, division, print_function

import os
from six.moves import cPickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import multigenAnalysisPlot

CRITICAL_N = [1, 2, 4, 8]


class Plot(multigenAnalysisPlot.MultigenAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, seedOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		sim_data = cPickle.load(open(simDataFile, "rb"))
		genomeLength = len(sim_data.process.replication.genome_sequence)

		ap = AnalysisPaths(seedOutDir, multi_gen_plot = True)

		# Get all cells
		allDir = ap.get_cells()

		# noinspection PyTypeChecker
		fig, axesList = plt.subplots(6, sharex = True)
		fig.set_size_inches(11, 11)
		for simDir in allDir:
			simOutDir = os.path.join(simDir, "simOut")
			time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") / 60. / 60.
			mass = TableReader(os.path.join(simOutDir, "Mass"))

			# Get fork positions
			replication_data_file = TableReader(os.path.join(simOutDir, "ReplicationData"))
			fork_coordinates = replication_data_file.readColumn("fork_coordinates")

			# Get fork counts
			pairsOfForks = np.logical_not(np.isnan(fork_coordinates)).sum(axis=1) / 2

			# Down sample dna polymerase position, every position is only plotted once here
			# using numpy ninja-ness
			unique, index, value = np.unique(fork_coordinates, return_index=True, return_inverse=True)
			m = np.zeros_like(value, dtype=bool)
			m[index] = True
			m = m.reshape(fork_coordinates.shape)
			fork_coordinates[~m] = np.nan

			# Skip plot if there are no replication forks in this generation
			if fork_coordinates.shape[1] > 0:
				axesList[0].plot(time, fork_coordinates, marker=',', markersize=1, linewidth=0)
			axesList[0].set_xlim([0, time.max()])
			axesList[0].set_yticks([-1 * genomeLength / 2, 0, genomeLength / 2])
			axesList[0].set_yticklabels(['-terC', 'oriC', '+terC'])
			axesList[0].set_ylabel("DNA polymerase\nposition (nt)")

			axesList[1].plot(time, pairsOfForks, linewidth=2)
			axesList[1].set_yticks(np.arange(0,7))
			axesList[1].set_ylim([0, 6])
			axesList[1].set_xlim([0, time.max()])
			axesList[1].plot([time.max(), time.max()], axesList[1].get_ylim(), 'k')
			axesList[1].set_ylabel("Pairs of\nforks")

			# Factors of critical initiation mass
			totalMass = mass.readColumn("cellMass")
			criticalInitiationMass = replication_data_file.readColumn("criticalInitiationMass")
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
			numberOfOric = replication_data_file.readColumn("numberOfOric")
			axesList[4].plot(time, numberOfOric, linewidth=2)
			axesList[4].set_ylabel("Number of\noriC")
			axesList[4].set_xlim([0, time.max()])

			# Mass per oriC
			criticalMassPerOriC = replication_data_file.readColumn("criticalMassPerOriC")
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
