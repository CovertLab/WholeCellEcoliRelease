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

import wholecell.utils.config

def cacheKnowledgeBase(outDir):
	sys.path.append(str(os.path.expanduser(wholecell.utils.config.KNOWLEDGEBASE_PACKAGE_DIR)))
	import ecoliwholecellkb_project.KnowledgeBaseEcoli

	kb = ecoliwholecellkb_project.KnowledgeBaseEcoli.KnowledgeBaseEcoli()
	cPickle.dump(kb, open(os.path.join(outDir, "KnowledgeBase.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

	return kb

def loadKnowledgeBase(inDir):
	sys.path.append(str(os.path.expanduser(wholecell.utils.config.KNOWLEDGEBASE_PACKAGE_DIR)))
	import ecoliwholecellkb_project.KnowledgeBaseEcoli
	kb = cPickle.load(open(inDir, "rb"))
	return kb
