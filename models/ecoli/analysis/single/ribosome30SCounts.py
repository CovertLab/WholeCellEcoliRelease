#!/usr/bin/env python
"""
Plots counts of 30S rRNA, associated proteins, and complexes

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/5/2014
"""

import argparse
import os

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

import wholecell.utils.constants

FONT = {'family' : 'normal',
		'weight':	'normal',
		'size'	:	8
		}

def setAxisMaxMin(axis, data):
	ymax = np.max(data)
	ymin = 0
	if ymin == ymax:
		axis.set_yticks([ymin])
	else:
		axis.set_yticks([ymin, ymax])

def sparklineAxis(axis, tickPos):
	axis.spines['top'].set_visible(False)
	axis.yaxis.set_ticks_position(tickPos)
	axis.xaxis.set_ticks_position('none')
	axis.tick_params(which = 'both', direction = 'out')
	axis.tick_params(labelbottom = 'off')


def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Load data from KB
	kb = cPickle.load(open(kbFile, "rb"))
	proteinIds = kb.s30_proteins
	rnaIds = [kb.monomerData['rnaId'][np.where(kb.monomerData['id'] == pid)[0][0]] for pid in proteinIds]
	rRnaIds = kb.s30_16sRRNA
	complexIds = [kb.s30_fullComplex]

	# Load count data for s30 proteins, rRNA, and final 30S complex
	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:
		# Get indexes
		moleculeIds = bulkMoleculesFile.root.names.moleculeIDs.read()
		proteinIndexes = np.array([moleculeIds.index(protein) for protein in proteinIds], np.int)
		rnaIndexes = np.array([moleculeIds.index(rna) for rna in rnaIds], np.int)
		rRnaIndexes = np.array([moleculeIds.index(rRna) for rRna in rRnaIds], np.int)
		complexIndexes = np.array([moleculeIds.index(comp) for comp in complexIds], np.int)

		# Load data
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules
		time = bulkMolecules.col("time")
		proteinCounts = bulkMolecules.read(0, None, 1, "counts")[:, proteinIndexes]
		rnaCounts = bulkMolecules.read(0, None, 1, "counts")[:, rnaIndexes]
		rRnaCounts = bulkMolecules.read(0, None, 1, "counts")[:, rRnaIndexes]
		complexCounts = bulkMolecules.read(0, None, 1, "counts")[:, complexIndexes]

	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	for idx in xrange(len(proteinIds)):
		rna_axis = plt.subplot(9, 3, idx + 1)

		rna_axis.step(time / 60., rnaCounts[:, idx], 'b', linewidth = 2)
		sparklineAxis(rna_axis, 'left')
		setAxisMaxMin(rna_axis, rnaCounts[:, idx])

		protein_axis = rna_axis.twinx()
		protein_axis.step(time / 60., proteinCounts[:, idx], 'r', linewidth = 2)
		sparklineAxis(protein_axis, 'right')		
		setAxisMaxMin(protein_axis, proteinCounts[:, idx])

		# Component label
		rna_axis.set_xlabel(proteinIds[idx][:-3])

	for idx in xrange(len(rRnaIds)):
		rna_axis = plt.subplot(9, 3, idx + len(proteinIds) + 1)

		rna_axis.step(time / 60., rRnaCounts[:, idx], 'b', linewidth = 2)
		sparklineAxis(rna_axis, 'left')
		setAxisMaxMin(rna_axis, rRnaCounts[:, idx])

		# Component label
		rna_axis.set_xlabel(rRnaIds[idx][:-3])

	for idx in xrange(len(complexIds)):
		complex_axis = plt.subplot(9, 3, idx + len(proteinIds) + len(rRnaIds) + 1)

		complex_axis.step(time / 60., complexCounts[:, idx], 'r', linewidth = 2)
		sparklineAxis(complex_axis, 'left')
		setAxisMaxMin(complex_axis, complexCounts[:, idx])

		# Component label
		complex_axis.set_xlabel(complexIds[idx][:-3])

	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))

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
