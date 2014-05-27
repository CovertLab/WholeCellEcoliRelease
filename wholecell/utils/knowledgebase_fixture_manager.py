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

_KB_FILENAME = "KnowledgeBase.cPickle"

def cacheKnowledgeBase(directory):
	kb = wholecell.reconstruction.knowledge_base_ecoli.KnowledgeBaseEcoli()
	
	writeKnowledgeBase(directory, kb)

	return kb


def loadKnowledgeBase(directory):
	kb = cPickle.load(open(os.path.join(directory, _KB_FILENAME), "rb"))
	return kb


def writeKnowledgeBase(directory, kb):
	cPickle.dump(kb, open(os.path.join(directory, _KB_FILENAME), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

