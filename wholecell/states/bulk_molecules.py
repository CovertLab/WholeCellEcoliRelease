#!/usr/bin/env python

"""
BulkMolecules.py

State which represents for a class of molecules the bulk copy numbers.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/04/2013
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

# TODO: break this into multiple files, it's becoming unbearably long

import re

import numpy as np
import tables

import wholecell.states.state as wcState
import wholecell.states.partition as wcPartition


class BulkMoleculesContainer(object):
	def __init__(self, moleculeNames):
		self._nMolecules = len(moleculeNames)

		self._moleculeNames = moleculeNames

		self._moleculeIndex = {moleculeName:index for index, moleculeName in enumerate(moleculeNames)}

		self._counts = np.zeros(self._nMolecules, np.uint64)


	def counts(self, names = None):
		if names == None:
			return self._counts

		else:
			return self._counts[self._namesToIndexes(names)]


	def countsIs(self, names, values):
		self._counts[self._namesToIndexes(names)] = values


	def countsInc(self, names, values): 
		self._counts[self._namesToIndexes(names)] += values


	def countsDec(self, names, values): # TODO: raise exception if > max?
		self._counts[self._namesToIndexes(names)] -= values


	def countsView(self, names):
		raise BulkMoleculesView(self, self._namesToIndexes(names))


	def count(self, name):
		return self._counts[self._moleculeIndex[name]]


	def countInc(self, name, value):
		self._counts[self._moleculeIndex[name]] += value


	def countDec(self, name, value): # TODO: raise exception if > max?
		self._counts[self._moleculeIndex[name]] -= value


	def countView(self, name):
		raise BulkMoleculeView(self, self._moleculeIndex[name])


	def _namesToIndexes(self, names):
		return np.array([self._moleculeIndex[name] for name in names])

	# TODO: mass calculation
	# TODO: saving


class BulkMoleculesView(object):
	def __init__(self, container, indexes):
		self._container = container
		self._indexes = indexes


	def counts(self):
		return self._container._counts[self._indexes]


	def countsIs(self, values):
		self._container._counts[self._indexes] = values


	def countsInc(self, values): 
		self._container._counts[self._indexes] += values


	def countsDec(self, values): # TODO: raise exception if > max?
		self._container._counts[self._indexes] -= values


class BulkMoleculeView(object):
	def __init__(self, container, index):
		self._container = container
		self._index = index


	def count(self):
		return self._container._counts[self._index]


	def countIs(self, values):
		self._container._counts[self._index] = values


	def countInc(self, values): 
		self._container._counts[self._index] += values


	def countDec(self, values): # TODO: raise exception if > max?
		self._container._counts[self._index] -= values



class BulkMolecules(wcState):
	pass


class BulkMoleculesPartition(wcPartition):
	pass
