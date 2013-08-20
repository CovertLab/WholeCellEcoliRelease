#!/usr/bin/env python

"""
generateTestFixtures
Generates fixtures for whole-cell tests

Example:
>>> from generateTestFixtures import *
>>> generateTestFixtures()

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/5/2013
"""

import os
import os.path
import cPickle
import pickle
import copy_reg
import types

# We need the following two methods so that we can pickle instancemethods
# Also need the copy_reg.pickle() line in generateTestFixtures()
# Borrowed from: http://mail.python.org/pipermail/python-list/2006-October/367078.html
# (Thanks Steven Bethard!)

def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

def generateTestFixtures():
	import wholecell.kb.KnowledgeBase
	import wholecell.sim.Simulation
	import wholecell.util.Fitter

	copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

	# Create output directory
	outDir = "data/fixtures"
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	# Construct KB
	kb = wholecell.kb.KnowledgeBase.KnowledgeBase(dataFileDir = "data/parsed/", seqFileName = "data/raw/sequence.txt")
	cPickle.dump(kb, open(os.path.join(outDir, "KnowledgeBase.cPickle"), "w"), protocol = cPickle.HIGHEST_PROTOCOL)

	# Construct simulation
	sim = wholecell.sim.Simulation.Simulation(kb)
	sim.setOptions({"seed": 1})
	wholecell.util.Fitter.Fitter.FitSimulation(sim, kb)
	sim.calcInitialConditions()
	cPickle.dump(sim, open(os.path.join(outDir, "Simulation.cPickle"), "w"), protocol = cPickle.HIGHEST_PROTOCOL)