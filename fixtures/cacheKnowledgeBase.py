#!/usr/bin/env python

"""
cacheKnowledgeBase
Generates cPickle file of knowledgebase for simulations

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/5/2014
"""

import os
import cPickle

import wholecell.reconstruction.knowledgebase

def main():
	# Create output directory
	outDir = "fixtures/sim"
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	# Construct KB
	kb = wholecell.reconstruction.knowledgebase.KnowledgeBase(dataFileDir = "data/parsed/", seqFileName = "data/raw/sequence.txt")
	cPickle.dump(kb, open(os.path.join(outDir, "KnowledgeBase.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
	main()
