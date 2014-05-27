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

FILENAME = "KnowledgeBase.cPickle"

def cacheKnowledgeBase(directory):
	kb = wholecell.reconstruction.knowledge_base_ecoli.KnowledgeBaseEcoli()
	cPickle.dump(kb, open(os.path.join(directory, FILENAME), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

	return kb

def loadKnowledgeBase(directory):
	kb = cPickle.load(open(os.path.join(directory, FILENAME), "rb"))
	return kb
