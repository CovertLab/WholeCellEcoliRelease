#!/usr/bin/env python

"""
createVariantKb.py

Modifies a knowledgebase.

"""

from __future__ import division

from models.ecoli.sim.variants.variants import nameToFunctionMapping

import argparse
import cPickle
import os

def main(variantFunction, index, inputKb, outputKb, metadataDir):

	if variantFunction not in nameToFunctionMapping:
		raise Exception, "%s is not a valid variantFunction!" % variantFunction

	kb = cPickle.load(open(inputKb, "rb"))

	info = nameToFunctionMapping[variantFunction](kb, index)

	print info["shortName"]

	cPickle.dump(
		kb,
		open(outputKb, "wb"),
		protocol = cPickle.HIGHEST_PROTOCOL
		)

	h = open(os.path.join(metadataDir, "short_name"), "w")
	h.write("%s\n" % info["shortName"])
	h.close()

	h = open(os.path.join(metadataDir, "description"), "w")
	h.write("%s\n" % info["desc"])
	h.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("variantFunction", help = "Function which modifies a kb object", type = str)
	parser.add_argument("index", help = "Index to pass to variant function", type = int)
	parser.add_argument("inputKb", help = "Input knowledgebase file", type = str)
	parser.add_argument("outputKb", help = "Output knowledgebase file", type = str)
	parser.add_argument("metadataDir", help = "Metadata directory", type = str)

	args = parser.parse_args().__dict__

	main(args["variantFunction"], args["index"], args["inputKb"], args["outputKb"], args["metadataDir"])