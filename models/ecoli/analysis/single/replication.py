#!/usr/bin/env python
"""
Plots things relevant to DNA replication

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/2015
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

PLACE_HOLDER = -1

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load KB
	kb = cPickle.load(open(kbFile, "rb"))

	oriC = kb.constants.oriCCenter.asNumber()
	terC = kb.constants.terCCenter.asNumber()
	genomeLength = len(kb.process.replication.genome_sequence)

	# Load time
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	# Load replication data
	dnaPolyFile = TableReader(os.path.join(simOutDir, "ReplicationData"))
	sequenceIdx = dnaPolyFile.readColumn("sequenceIdx")
	sequenceLength = dnaPolyFile.readColumn("sequenceLength")
	numberOfOric = dnaPolyFile.readColumn("numberOfOric")
	criticalMassPerOriC = dnaPolyFile.readColumn("criticalMassPerOriC")
	dnaPolyFile.close()

	# Load dna mass data
	massFile = TableReader(os.path.join(simOutDir, "Mass"))
	totalMass = massFile.readColumn("cellMass")
	dnaMass = massFile.readColumn("dnaMass")
	massFile.close()

	# Setup elongation length data
	reverseIdx = 1
	reverseCompIdx = 3
	reverseSequences = np.logical_or(sequenceIdx == reverseIdx, sequenceIdx == reverseCompIdx)
	sequenceLength[reverseSequences] = -1 * sequenceLength[reverseSequences]
	sequenceLength[sequenceLength == PLACE_HOLDER] = np.nan

	# Count pairs of forks, initation, and termination events
	pairsOfForks = (sequenceIdx != PLACE_HOLDER).sum(axis = 1) / 4

	# Count chromosome equivalents
	chromMass = (kb.getter.getMass(['CHROM_FULL[c]'])[0] / kb.constants.nAvogadro).asNumber(units.fg)
	chromEquivalents = dnaMass / chromMass

	# Count full chromosomes
	bulkMoleculesFile = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	bulkIds = bulkMoleculesFile.readAttribute("objectNames")
	chromIdx = bulkIds.index("CHROM_FULL[c]")
	fullChromosomeCounts = bulkMoleculesFile.readColumn("counts")[:,chromIdx]

	# Count 60 min doubling time mass equivalents
	avgCell60MinDoublingTimeTotalMassInit = kb.mass.avgCell60MinDoublingTimeTotalMassInit.asNumber(units.fg)
	sixtyMinDoublingInitMassEquivalents = totalMass / avgCell60MinDoublingTimeTotalMassInit

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
	ax.plot(time / 60., sixtyMinDoublingInitMassEquivalents, linewidth=2)
	ax.set_xticks([0, time.max() / 60])
	ax.set_yticks(np.arange(1., 8., 0.5))
	ax.set_ylim([np.around(sixtyMinDoublingInitMassEquivalents.min(), decimals=1) - 0.1, np.around(sixtyMinDoublingInitMassEquivalents.max(), decimals=1) + 0.1])
	ax.set_ylabel("Equivalents of initial\nmass for $t_d=60$ min")
	
	ax = plt.subplot(7,1,5, sharex=ax)
	ax.plot(time / 60., criticalMassPerOriC, linewidth=2)
	ax.set_xticks([0, time.max() / 60])
	ax.set_yticks([0.5, 1.0])
	# ax.set_ylim([0.4, 1.1])
	ax.set_ylabel("Critical mass\nper oriC")

	ax = plt.subplot(7,1,6, sharex=ax)
	ax.plot(time / 60., numberOfOric, linewidth=2)
	ax.set_xticks([0, time.max() / 60])
	ax.set_ylabel("Number of\noriC")
	ax.set_ylim([0, numberOfOric.max() + 1])

	ax = plt.subplot(7,1,7, sharex=ax)
	ax.plot(time / 60., fullChromosomeCounts, linewidth=2)
	ax.set_xticks([0, time.max() / 60])
	ax.set_ylabel("Full\nchromosomes")
	ax.set_ylim([0, fullChromosomeCounts.max() + 1])	

	ax.set_xlim([0, time.max() / 60])
	ax.set_xlabel("Time (min)")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")

	plt.figure(figsize = (8.5, 11))

	for idx in [0,1,2,3]:
		ax = plt.subplot(4,1,idx+1)
		data = (sequenceIdx == idx).sum(axis=1)
		ax.plot(time / 60., data)
		ax.set_ylim([0, data.max()+1])

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, 'replicationSequenceIdx')
	plt.close("all")

if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])
