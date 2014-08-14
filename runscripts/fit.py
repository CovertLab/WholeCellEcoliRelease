#!/usr/bin/env python

"""
fit.py

Fit knowledgebase fixture.
"""

import cPickle
import os
import argparse

import reconstruction.ecoli.knowledge_base
import reconstruction.ecoli.fitter
import wholecell.utils.constants

def main(fitLevel, kbDirectory):
	kbDirectory = kbDirectory

	prevFitLevel = fitLevel - 1

	prevFitKbFileName = (
		wholecell.utils.constants.SERIALIZED_KB_PREFIX +
		("_Fit_%d" % prevFitLevel) +
		wholecell.utils.constants.SERIALIZED_KB_SUFFIX
		)

	kb = cPickle.load(open(os.path.join(kbDirectory, prevFitKbFileName), "rb"))

	simOutDir = None

	print "Fitting knowledgebase (Level %d)" % fitLevel

	reconstruction.ecoli.fitter.fitAtLevel(fitLevel, kb, simOutDir)

	fileName = (
				wholecell.utils.constants.SERIALIZED_KB_PREFIX +
				("_Fit_%d" % fitLevel) +
				wholecell.utils.constants.SERIALIZED_KB_SUFFIX
				)

	print "Saving fit knowledgebase"
	cPickle.dump(
		kb,
		open(os.path.join(
			kbDirectory,
			fileName,
			), "wb"),
		protocol = cPickle.HIGHEST_PROTOCOL
		)

	# Create symlink indicating that this is the most fit knowledgebase
	# Note: this assumes things are getting called in order, etc.
	symlink = os.path.join(
		kbDirectory,
		wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME		
		)

	if os.path.exists(symlink):
		os.unlink(symlink)

	os.symlink(fileName, symlink)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"fitLevel",
		help = "Stage of fitting to perform",
		type = int)
	parser.add_argument(
		"kbDirectory",
		help = "Directory containing output",
		type = str)

	args = parser.parse_args().__dict__

	main(args["fitLevel"], args["kbDirectory"])
