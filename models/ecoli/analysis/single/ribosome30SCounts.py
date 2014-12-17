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
from wholecell.utils.sparkline import sparklineAxis, setAxisMaxMinY

FONT = {
		'size'	:	8
		}

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
		freeProteinCounts = bulkMolecules.read(0, None, 1, "counts")[:, proteinIndexes]
		rnaCounts = bulkMolecules.read(0, None, 1, "counts")[:, rnaIndexes]
		freeRRnaCounts = bulkMolecules.read(0, None, 1, "counts")[:, rRnaIndexes]
		complexCounts = bulkMolecules.read(0, None, 1, "counts")[:, complexIndexes]

	with tables.open_file(os.path.join(simOutDir, "UniqueMoleculeCounts.hdf")) as uniqueMoleculesFile:
		uniqueMoleculeCounts = uniqueMoleculesFile.root.UniqueMoleculeCounts
		ribosomeIndex = uniqueMoleculeCounts.attrs.uniqueMoleculeIds.index("activeRibosome")
		activeRibosome = uniqueMoleculeCounts.col("uniqueMoleculeCounts")[:, ribosomeIndex]

	plt.figure(figsize = (8.5, 11))
	matplotlib.rc('font', **FONT)

	for idx in xrange(len(proteinIds)):
		rna_axis = plt.subplot(12, 3, idx + 1)

		sparklineAxis(rna_axis, time / 60., rnaCounts[:, idx], 'left', '-', 'b')
		setAxisMaxMinY(rna_axis, rnaCounts[:, idx])

		protein_axis = rna_axis.twinx()
		sparklineAxis(protein_axis, time / 60., freeProteinCounts[:, idx], 'right', '-', 'r')
		setAxisMaxMinY(protein_axis, freeProteinCounts[:, idx])

		# Component label
		rna_axis.set_xlabel(proteinIds[idx][:-3])

	for idx in xrange(len(rRnaIds)):
		rna_axis = plt.subplot(12, 3, idx + len(proteinIds) + 1)

		sparklineAxis(rna_axis, time / 60., freeRRnaCounts[:, idx], 'left', '-', 'b')
		setAxisMaxMinY(rna_axis, freeRRnaCounts[:, idx])

		# Component label
		rna_axis.set_xlabel(rRnaIds[idx][:-3])

	for idx in xrange(len(complexIds)):
		complex_axis = plt.subplot(12, 3, idx + len(proteinIds) + len(rRnaIds) + 1)

		sparklineAxis(complex_axis, time / 60., complexCounts[:, idx], 'left', '-', 'r')
		setAxisMaxMinY(complex_axis, complexCounts[:, idx])

		# Component label
		complex_axis.set_xlabel(complexIds[idx][:-3])

	# Plot number of ribosomes
	ribosome_axis = plt.subplot(12, 3, 1 + len(proteinIds) + len(rRnaIds) + len(complexIds) + 1)
	sparklineAxis(ribosome_axis, time / 60., activeRibosome, 'left', '-', 'r')
	setAxisMaxMinY(ribosome_axis, activeRibosome)
	ribosome_axis.set_xlabel('Active ribosome')

	# Save
	plt.subplots_adjust(hspace = 0.5, wspace = 0.5)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)

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
