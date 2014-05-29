#!/usr/bin/env python

"""
createKbs.py

Create unfit and fit knowledgebase fixtures.
"""

import wholecell.reconstruction.knowledge_base_ecoli
import wholecell.reconstruction.fitter
import wholecell.utils.config
import cPickle
import os
import argparse

_UNFIT_FILENAME = "KnowledgeBase_Unfit.cPickle"
_FIT_FILENAME = "KnowledgeBase_Fit.cPickle"

def main(outputDirectory):

	if not os.path.exists(outputDirectory):
		os.makedirs(outputDirectory)

	print "Instantiating unfit knowledgebase"
	kb = wholecell.reconstruction.knowledge_base_ecoli.KnowledgeBaseEcoli()

	print "Saving unfit knowledgebase"
	cPickle.dump(
		kb,
		open(os.path.join(
			outputDirectory,
			wholecell.utils.config.SERIALIZED_KB_UNFIT_FILENAME
			), "wb"),
		protocol = cPickle.HIGHEST_PROTOCOL
		)

	print "Fitting knowledgebase"
	wholecell.reconstruction.fitter.fitKb(kb)

	print "Saving fit knowledgebase"
	cPickle.dump(
		kb,
		open(os.path.join(
			outputDirectory,
			wholecell.utils.config.SERIALIZED_KB_FIT_FILENAME
			), "wb"),
		protocol = cPickle.HIGHEST_PROTOCOL
		)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("outputDirectory", help = "Directory containing output", type = str)

	args = parser.parse_args().__dict__

	main(args["outputDirectory"])