#!/usr/bin/env python

"""
getNumVariants.py

Get the number of variants for a specific variant function.

"""

from __future__ import division

from models.ecoli.sim.variants.variants import nameToNumIndicesMapping

import argparse
import cPickle

def main(variantFunction, inputKb):

	if variantFunction not in nameToNumIndicesMapping:
		raise Exception, "%s is not a valid variantFunction!" % variantFunction

	kb = cPickle.load(open(inputKb, "rb"))

	print "%d" % (nameToNumIndicesMapping[variantFunction](kb))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("variantFunction", help = "Function which modifies a kb object", type = str)
	parser.add_argument("inputKb", help = "Input knowledgebase file", type = str)

	args = parser.parse_args().__dict__

	main(args["variantFunction"], args["inputKb"])