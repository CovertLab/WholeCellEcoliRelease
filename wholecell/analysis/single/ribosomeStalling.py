#!/usr/bin/env python
"""
Plot ribosome stalling

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/22/2014
"""

from __future__ import division

import argparse
import os

import tables
import numpy as np
import matplotlib
# matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.utils.knowledgebase_fixture_manager import loadKnowledgeBase
from wholecell.utils.constants import SIM_FIXTURE_DIR

def main(simOutDir, plotOutDir, plotOutFileName):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	with tables.open_file(os.path.join(simOutDir, "RibosomeStalling.hdf")) as h5file:
		ribosomeStalling = h5file.root.RibosomeStalling.read()

	times, inverseTimeIndexes = np.unique(ribosomeStalling["t"], return_inverse = True)

	stallsByTime = []

	for i, time in enumerate(times):
		stallsByTime.append(
			ribosomeStalling["stalls"][np.where(i == inverseTimeIndexes)]
			)

	proteins, inverseProteinIndexes = np.unique(ribosomeStalling["proteinIndex"], return_inverse = True)

	stallsByProtein = []

	for i, protein in enumerate(proteins):
		stallsByProtein.append(
			ribosomeStalling["stalls"][np.where(i == inverseProteinIndexes)]
			)

	kb = loadKnowledgeBase(os.path.join(SIM_FIXTURE_DIR, "KnowledgeBase.cPickle"))

	proteinNames = kb.monomerData["id"]
	aaCounts = kb.monomerData["aaCounts"].magnitude

	aaCounts = np.delete(aaCounts, kb.aaIDs.index("SEC-L[c]"), 1)

	proteinComposition = aaCounts / aaCounts.sum(1)[:, np.newaxis]

	metaboliteIds = kb.wildtypeBiomass["metaboliteId"].tolist()

	aaCountsBiomass = kb.wildtypeBiomass["biomassFlux"][np.array([
		metaboliteIds.index(aa) for aa in kb.aaIDs if aa in metaboliteIds
		])].magnitude

	averageProteinComposition = aaCountsBiomass / aaCountsBiomass.sum()

	compositionDeviation = np.sqrt(
		((proteinComposition - averageProteinComposition)**2).sum(1)
		)[proteins]

	totalStallsByTime = np.array([
		values.sum() for values in stallsByTime
		])

	meanStallsByTime = np.array([
		values.mean() for values in stallsByTime
		])

	meanStallsByStalledByTime = np.array([
		values[values > 0].mean() for values in stallsByTime
		])

	stdStallsByTime = np.array([
		values.std() for values in stallsByTime
		])

	fractionStalledByTime = np.array([
		(values > 0).sum() / values.size for values in stallsByTime
		])

	totalStallsByProtein = np.array([
		values.sum() for values in stallsByProtein
		])

	meanStallsByProtein = np.array([
		values.mean() for values in stallsByProtein
		])

	meanStallsByStalledByProtein = np.array([
		values[values > 0].mean() for values in stallsByProtein
		])

	stdStallsByProtein = np.array([
		values.std() for values in stallsByProtein
		])

	fractionStalledByProtein = np.array([
		(values > 0).sum() / values.size for values in stallsByProtein
		])

	# Plots by time

	# plt.figure(figsize = (8.5, 11))

	plt.subplot(3, 1, 0)

	plt.plot(times / 60, totalStallsByTime)
	plt.title("Total number of stalls, by time")

	plt.xlabel("Time (min)")
	plt.ylabel("Stalls")

	plt.subplot(3, 1, 1)

	plt.plot(times / 60, meanStallsByStalledByTime)
	# plt.plot(times / 60, meanStallsByTime + stdStallsByTime, '--')
	# plt.plot(times / 60, meanStallsByTime - stdStallsByTime, '--')

	plt.title("Average number of stalls, by time")

	plt.xlabel("Time (min)")
	plt.ylabel("Stalls")

	plt.subplot(3, 1, 2)

	plt.plot(times / 60, fractionStalledByTime)
	plt.title("Fraction stalled, by time")

	plt.xlabel("Time (min)")
	plt.ylabel("Stalls")

	# plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	plt.show()

	# Plots by protein

	# plt.figure(figsize = (8.5, 11))

	plt.subplot(3, 1, 0)

	plt.plot(proteins, totalStallsByProtein)
	plt.title("Total number of stalls, by protein")

	plt.xlabel("Protein index")
	plt.ylabel("Stalls")

	plt.subplot(3, 1, 1)

	plt.plot(proteins, meanStallsByStalledByProtein)
	# plt.plot(proteins, meanStallsByProtein + stdStallsByProtein, '--')
	# plt.plot(proteins, meanStallsByProtein - stdStallsByProtein, '--')

	plt.title("Average number of stalls, by protein")

	plt.xlabel("Protein index")
	plt.ylabel("Stalls")

	plt.subplot(3, 1, 2)

	plt.plot(proteins, fractionStalledByProtein)
	plt.title("Fraction stalled, by protein")

	plt.xlabel("Protein index")
	plt.ylabel("Stalls")

	# # plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	plt.show()

	# plt.savefig(os.path.join(plotOutDir, plotOutFileName))

	import ipdb; ipdb.set_trace()


if __name__ == "__main__":
	# parser = argparse.ArgumentParser()
	# parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	# parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	# parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	# args = parser.parse_args().__dict__

	# main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])

	# TODO: actually use this code/output something helpful
	pass
