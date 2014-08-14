#!/usr/bin/env python

"""
buildKb.py

Create unfit knowledgebase fixture.
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
	fileName = (
				wholecell.utils.constants.SERIALIZED_KB_PREFIX +
				"_Fit_0" +
				wholecell.utils.constants.SERIALIZED_KB_SUFFIX
				)

	cPickle.dump(
		kb,
		open(os.path.join(
			outputDirectory,
			fileName
			), "wb"),
		protocol = cPickle.HIGHEST_PROTOCOL
		)

	# Create symlink indicating that this is the unfit knowledgebase
	symlink = os.path.join(
		outputDirectory,
		wholecell.utils.constants.SERIALIZED_KB_UNFIT_FILENAME
		)

	if os.path.lexists(symlink):
		os.unlink(symlink)

	os.symlink(fileName, symlink)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"--outputDirectory",
		help = "Directory containing output",
		type = str)

	args = parser.parse_args().__dict__

	main(args["outputDirectory"])
