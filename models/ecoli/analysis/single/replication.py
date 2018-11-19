"""
Plots simulation outputs relevant to DNA replication

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/2015
"""

from __future__ import absolute_import, division, print_function

import os

import numpy as np
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
from wholecell.utils import units
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot

PLACE_HOLDER = -1

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(simOutDir):
			raise Exception, "simOutDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		# Load KB
		sim_data = cPickle.load(open(simDataFile, "rb"))

		oriC = sim_data.constants.oriCCenter.asNumber()
		terC = sim_data.constants.terCCenter.asNumber()
		genomeLength = len(sim_data.process.replication.genome_sequence)

		# Load time
		initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
		time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

		# Load replication data
		dnaPolyFile = TableReader(os.path.join(simOutDir, "ReplicationData"))
		sequenceIdx = dnaPolyFile.readColumn("sequenceIdx")
		sequenceLength = dnaPolyFile.readColumn("sequenceLength")
		numberOfOric = dnaPolyFile.readColumn("numberOfOric")
		criticalMassPerOriC = dnaPolyFile.readColumn("criticalMassPerOriC")
		criticalInitiationMass = dnaPolyFile.readColumn("criticalInitiationMass")

		# Load dna mass data
		massFile = TableReader(os.path.join(simOutDir, "Mass"))
		totalMass = massFile.readColumn("cellMass")
		dnaMass = massFile.readColumn("dnaMass")

		# Setup elongation length data
		reverseIdx = 1
		reverseCompIdx = 3
		reverseSequences = np.logical_or(sequenceIdx == reverseIdx, sequenceIdx == reverseCompIdx)
		sequenceLength[reverseSequences] = -1 * sequenceLength[reverseSequences]
		sequenceLength[sequenceLength == PLACE_HOLDER] = np.nan

		# Count pairs of forks, initiation, and termination events
		pairsOfForks = (sequenceIdx != PLACE_HOLDER).sum(axis = 1) / 4

		# Count chromosome equivalents
		chromMass = (sim_data.getter.getMass(['CHROM_FULL[c]'])[0] / sim_data.constants.nAvogadro).asNumber(units.fg)
		chromEquivalents = dnaMass / chromMass

		# Count full chromosomes
		unique_molecule_counts_reader = TableReader(
			os.path.join(simOutDir, "UniqueMoleculeCounts"))
		full_chromosome_index = unique_molecule_counts_reader.readAttribute(
			"uniqueMoleculeIds").index("fullChromosome")
		full_chromosome_counts = unique_molecule_counts_reader.readColumn(
			"uniqueMoleculeCounts")[:, full_chromosome_index]

		# Count critical initiation mass equivalents
		# criticalInitiationMass[0] = criticalInitiationMass[1]
		criticalMassEquivalents = totalMass / criticalInitiationMass

		# Plot stuff
		plt.figure(figsize = (8.5, 11))

		ax = plt.subplot(7,1,1)
		ax.plot(time / 60., sequenceLength, marker='.', markersize=1, linewidth=0)
		ax.set_xticks([0, time.max() / 60])
		ax.set_yticks([-1 * genomeLength / 2, 0, genomeLength / 2])
		ax.set_yticklabels(['-terC', 'oriC', '+terC'])
		ax.set_ylabel("DNA polymerase\nposition (nt)")

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
