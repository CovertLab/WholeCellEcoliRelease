#!/usr/bin/env python

"""
knowledgebase_fixture_manager.py
Generates cPickle file of knowledgebase for simulations and tests

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/5/2014
"""

import os
import sys
import cPickle

import wholecell.reconstruction.knowledge_base_ecoli

def cacheKnowledgeBase(outDir):
	kb = wholecell.reconstruction.knowledge_base_ecoli.KnowledgeBaseEcoli()
	cPickle.dump(kb, open(os.path.join(outDir, "KnowledgeBase.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

	return kb

def loadKnowledgeBase(inDir):
	kb = cPickle.load(open(inDir, "rb"))
	return kb
