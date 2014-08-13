#!/usr/bin/env python

"""
createKbs.py

Create unfit and fit knowledgebase fixtures.
"""

import cPickle
import os
import argparse

import reconstruction.ecoli.knowledge_base
import reconstruction.ecoli.fitter
import wholecell.utils.constants

def main(outputDirectory = None):
	outputDirectory = outputDirectory or wholecell.utils.constants.SERIALIZED_KB_DIR

	if not os.path.exists(outputDirectory):
		os.makedirs(outputDirectory)

	print "Instantiating unfit knowledgebase"
	kb = reconstruction.ecoli.knowledge_base.KnowledgeBaseEcoli()

	print "Saving unfit knowledgebase"
	cPickle.dump(
		kb,
		open(os.path.join(
			outputDirectory,
			wholecell.utils.constants.SERIALIZED_KB_UNFIT_FILENAME
			), "wb"),
		protocol = cPickle.HIGHEST_PROTOCOL
		)

	print "Fitting knowledgebase"
	reconstruction.ecoli.fitter.fitKb(kb)

	print "Saving fit knowledgebase"
	cPickle.dump(
		kb,
		open(os.path.join(
			outputDirectory,
			wholecell.utils.constants.SERIALIZED_KB_FIT_FILENAME
			), "wb"),
		protocol = cPickle.HIGHEST_PROTOCOL
		)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"--outputDirectory",
		help = "Directory containing output",
		type = str)

	args = parser.parse_args().__dict__

	main(args["outputDirectory"])
