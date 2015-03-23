#!/usr/bin/env python

"""
buildKb.py

Create unfit knowledgebase fixture.
"""

import cPickle
import os
import argparse

import reconstruction.ecoli.simulation_data
import wholecell.utils.constants

def main(kbDirectory):

	if not os.path.exists(kbDirectory):
		os.makedirs(kbDirectory)

	print "Instantiating unfit knowledgebase"
	kb = reconstruction.ecoli.simulation_data.SimulationDataEcoli()

	print "Saving unfit knowledgebase"
	fileName = (
				wholecell.utils.constants.SERIALIZED_KB_PREFIX +
				"_Fit_0" +
				wholecell.utils.constants.SERIALIZED_KB_SUFFIX
				)

	cPickle.dump(
		kb,
		open(os.path.join(
			kbDirectory,
			fileName
			), "wb"),
		protocol = cPickle.HIGHEST_PROTOCOL
		)

	# Create symlink indicating that this is the unfit knowledgebase
	symlink = os.path.join(
		kbDirectory,
		wholecell.utils.constants.SERIALIZED_KB_UNFIT_FILENAME
		)

	if os.path.exists(symlink):
		os.unlink(symlink)

	os.symlink(fileName, symlink)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"kbDirectory",
		help = "Directory containing output kb",
		type = str)

	args = parser.parse_args().__dict__

	main(args["kbDirectory"])
