"""
Plots simulation outputs relevant to DNA replication
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
from six.moves import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		# Load KB
		sim_data = cPickle.load(open(simDataFile, "rb"))

		genomeLength = len(sim_data.process.replication.genome_sequence)

		# Load time
		main_reader = TableReader(os.path.join(simOutDir, "Main"))
		initialTime = main_reader.readAttribute("initialTime")
		time = main_reader.readColumn("time") - initialTime

		# Load replication data
		replication_data_file = TableReader(os.path.join(simOutDir, "ReplicationData"))
		fork_coordinates = replication_data_file.readColumn("fork_coordinates")
		numberOfOric = replication_data_file.readColumn("numberOfOric")
		criticalMassPerOriC = replication_data_file.readColumn("criticalMassPerOriC")
		criticalInitiationMass = replication_data_file.readColumn("criticalInitiationMass")

		# Load dna mass data
		massFile = TableReader(os.path.join(simOutDir, "Mass"))
		totalMass = massFile.readColumn("cellMass")
		dnaMass = massFile.readColumn("dnaMass")

		# Count pairs of forks, initiation, and termination events
		pairsOfForks = np.logical_not(np.isnan(fork_coordinates)).sum(axis = 1)/2

		# Count chromosome equivalents
		chromMass = (sim_data.getter.get_mass([sim_data.molecule_ids.full_chromosome])[0] / sim_data.constants.n_avogadro).asNumber(units.fg)
		chromEquivalents = dnaMass / chromMass

		# Count full chromosomes
		unique_molecule_counts_reader = TableReader(
			os.path.join(simOutDir, "UniqueMoleculeCounts"))
		full_chromosome_index = unique_molecule_counts_reader.readAttribute(
			"uniqueMoleculeIds").index('full_chromosome')
		full_chromosome_counts = unique_molecule_counts_reader.readColumn(
			"uniqueMoleculeCounts")[:, full_chromosome_index]

		# Count critical initiation mass equivalents
		if np.all(criticalInitiationMass > 0):
			criticalMassEquivalents = totalMass / criticalInitiationMass
		else:  # Cell does not have any chromosomes
			criticalMassEquivalents = np.zeros_like(totalMass)

		# Plot stuff
		plt.figure(figsize = (8.5, 11))

		ax = plt.subplot(7,1,1)
		# Skip if there are no replication forks
		if fork_coordinates.shape[1] > 0:
			ax.plot(time / 60., fork_coordinates, marker='.', markersize=1, linewidth=0)
		ax.set_xticks([0, time.max() / 60])
		ax.set_yticks([-1 * genomeLength / 2, 0, genomeLength / 2])
		ax.set_yticklabels(['-terC', 'oriC', '+terC'])
		ax.set_ylabel("Fork\nposition (nt)")

		ax = plt.subplot(7,1,2, sharex=ax)
		ax.plot(time / 60., chromEquivalents, linewidth=2)
		ax.set_xticks([0, time.max() / 60])
		ax.set_yticks(np.arange(chromEquivalents.min(), chromEquivalents.max() + 0.5, 0.5))
		ax.set_ylabel("Chromosome\nequivalents")

		ax = plt.subplot(7,1,3, sharex=ax)
		ax.plot(time / 60., pairsOfForks, linewidth=2)
		ax.set_xticks([0, time.max() / 60])
		ax.set_yticks(np.arange(0,7))
		ax.set_ylim([0, 6])
		ax.set_ylabel("Pairs of forks")

		ax = plt.subplot(7,1,4, sharex=ax)
		ax.plot(time / 60., criticalMassEquivalents, linewidth=2)
		ax.set_xticks([0, time.max() / 60])
		ax.set_yticks(np.arange(1., 8., 0.5))
		ax.set_ylim([np.around(criticalMassEquivalents[1:].min(), decimals=1) - 0.1, np.around(criticalMassEquivalents[1:].max(), decimals=1) + 0.1])
		ax.set_ylabel("Factors of critical\ninitiation mass")

		ax = plt.subplot(7,1,5, sharex=ax)
		ax.plot(time / 60., criticalMassPerOriC, linewidth=2)
		ax.plot(time / 60., np.ones_like(time), "k--", linewidth=2)
		ax.set_xticks([0, time.max() / 60])
		ax.set_yticks([0.5, 1.0])
		ax.set_ylabel("Critical mass\nper oriC")

		ax = plt.subplot(7,1,6, sharex=ax)
		ax.plot(time / 60., numberOfOric, linewidth=2)
		ax.set_xticks([0, time.max() / 60])
		ax.set_ylabel("Number of\noriC")
		ax.set_ylim([0, numberOfOric.max() + 1])

		ax = plt.subplot(7,1,7, sharex=ax)
		ax.plot(time / 60., full_chromosome_counts, linewidth=2)
		ax.set_xticks([0, time.max() / 60])
		ax.set_ylabel("Full\nchromosomes")
		ax.set_ylim([0, full_chromosome_counts.max() + 1])

		ax.set_xlim([0, time.max() / 60])
		ax.set_xlabel("Time (min)")

		plt.tight_layout()

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")


if __name__ == "__main__":
	Plot().cli()
