#!/usr/bin/env python

"""
uniques.py

State which represents for a class of molecules the bulk copy numbers as an array
and unique instances in a series of lists sharing a common index.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/04/2013
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

import numpy
# import types
# import sys
# import time

# TODO: Change 'object' to wholecell.sim.state.State
class MoleculesContainer(object):

	def __init__(self):
		self._countsBulk = None
		self._countsUnique = None
		self._massSingle = None
		self._dmass = None
		self._uniqueDict = None

		self._widIdx = None
		self._cIdx = None

		self._molecules = {}

	def initialize(self, kb):
		self._countsBulk = 0.0 * numpy.ones((len(kb.molecules), len(kb.compartments)))
		self._countsUnique = 0.0 * numpy.ones((len(kb.molecules), len(kb.compartments)))
		self._massSingle = numpy.array([[x["mass"]] * len(kb.compartments) for x in kb.molecules])

		self._uniqueDict = []
		for mol in kb.molecules:
			if mol["uniqueAttrs"] != None:
				self._uniqueDict.append([dict(dict(zip(mol["uniqueAttrs"] + ["objects"], [[] for x in xrange(len(mol["uniqueAttrs"]) + 1)]))) for x in kb.compartments])
			else:
				self._uniqueDict.append([{} for x in kb.compartments])

		self._dmass = numpy.zeros((len(kb.molecules), len(kb.compartments)))

		self._widIdx = dict([(x[0]["id"], x[1]) for x in zip(kb.molecules, range(len(kb.molecules)))])
		self._cIdx = dict([(x[0]["id"], x[1]) for x in zip(kb.compartments, range(len(kb.compartments)))])

	def molecule(self, wid, comp):
		if (wid, comp) not in self._molecules:
			self._molecules[wid, comp] = _Molecule(self, self._widIdx[wid], self._cIdx[comp], wid)

		return self._molecules[wid, comp]

def _uniqueInit(self, uniqueIdx):
	self._uniqueIdx = uniqueIdx

def _makeGetter(attr):
	def attrGetter(self):
		molCont, uniqueIdx = self._container, self._uniqueIdx
		molRowIdx, molColIdx = self._molRowIdx, self._molColIdx

		return molCont._uniqueDict[molRowIdx][molColIdx][attr][uniqueIdx]

	return attrGetter

def _makeSetter(attr):
	def attrSetter(self, newVal):
		molCont, uniqueIdx = self._container, self._uniqueIdx
		molRowIdx, molColIdx = self._molRowIdx, self._molColIdx

		molCont._uniqueDict[molRowIdx][molColIdx][attr][uniqueIdx] = newVal

	return attrSetter

class _Molecule(object):
	uniqueClassRegistry = {}
	def __init__(self, container, rowIdx, colIdx, wid):
		self._container = container
		self._rowIdx = rowIdx
		self._colIdx = colIdx
		self._wid = wid

		if len(self._container._uniqueDict[self._rowIdx][self._colIdx]) == 0:
			# Molecule has no attributes
			pass

		else:
			if self._wid in self.uniqueClassRegistry:
				# Molecule has been registered as a unique instance; use that class definition
				self._MoleculeUnique = self.uniqueClassRegistry[self._wid]
				self._MoleculeUnique._container = self._container
				self._MoleculeUnique._molRowIdx = self._rowIdx
				self._MoleculeUnique._molColIdx = self._colIdx

				for attr in self._container._uniqueDict[self._rowIdx][self._colIdx]:
					if not hasattr(self._MoleculeUnique, attr):
						setattr(self._MoleculeUnique, attr, _makeGetter(attr))

					if not hasattr(self._MoleculeUnique, attr + "Is"):
						setattr(self._MoleculeUnique, attr + "Is", _makeGetter(attr + "Is"))

			else:
				# Molecule has unique attributes, but isn't registered
				uniqueClassDefDict = {}
				uniqueClassDefDict["_container"] = self._container
				uniqueClassDefDict["_molRowIdx"] = self._rowIdx
				uniqueClassDefDict["_molColIdx"] = self._colIdx

				uniqueClassDefDict["__init__"] = _uniqueInit

				for attr in self._container._uniqueDict[self._rowIdx][self._colIdx]:
					uniqueClassDefDict[attr] = _makeGetter(attr)
					uniqueClassDefDict[attr + "Is"] = _makeSetter(attr)

				self._MoleculeUnique = type("MoleculeUnique", (), uniqueClassDefDict)

	def countsBulk(self):
		# Returns bulk count of molecule as a float
		return self._container._countsBulk[self._rowIdx, self._colIdx]

	def countsBulkIs(self, newVal):
		# Sets bulk count of molecule
		self._container._countsBulk[self._rowIdx, self._colIdx] = newVal

	def countsBulkInc(self, incVal):
		# Increments counts bulk by incVal
		self._container._countsBulk[self._rowIdx, self._colIdx] += incVal

	def countsBulkDec(self, decVal):
		# Decrements counts bulk by decVal
		self.countsBulkInc(-1 * decVal)

	def countsUnique(self):
		# Returns unique count of molecule as a float
		return self._container._countsUnique[self._rowIdx, self._colIdx]
	
	def dMassIs(self, newVal):
		# Sets the value for dMass
		self._container._dmass[self._rowIdx, self._colIdx] = newVal

	def dMassInc(self, incVal):
		# Increments dmass by incVal
		self._container._dmass[self._rowIdx, self._colIdx] += incVal

	def dMassDec(self, decVal):
		# Decrements count of dmass by decVal
		self.dMassInc(-1 * decVal)

	def massSingle(self):
		# Returns mass of single object
		return self._container._massSingle[self._rowIdx, self._colIdx]

	def massAll(self):
		# Returns mass of all objects bulk and unique
		return (self.countsBulk() + self.countsUnique()) * self.massSingle() + self._container._dmass[self._rowIdx, self._colIdx]

	def uniqueNew(self, attrs = None):
		# Creates new unique object with attributes defined by attrs
		# attrs should be in format: {"attr1" : value1, "attr2" : value2, ...}
		uniqueDict = self._container._uniqueDict[self._rowIdx][self._colIdx]
		
		if not len(uniqueDict):
			raise uniqueException, 'Attempting to create unique from object with no unique attributes!\n'
		if attrs != None and len(set(attrs).difference(set(uniqueDict.keys()))):
			raise uniqueException, 'A specified attribute is not included in knoweldge base for this unique object!\n'

		for attr in uniqueDict:
			if attrs != None and attr in attrs:
				uniqueDict[attr].append(attrs[attr])
			else:
				uniqueDict[attr].append(None)
		uniqueIdx = len(uniqueDict["objects"]) - 1
		uniqueDict["objects"][uniqueIdx] = self._MoleculeUnique(uniqueIdx)
		self._container._countsUnique[self._rowIdx, self._colIdx] += 1
		return uniqueDict["objects"][uniqueIdx]

	def uniquesWithAttrs(self, attrs = None):
		# Returns list of objects with attributes specified in attrs
		# attrs should be in format: {"attr1" : value1, "attr2" : value2, ...}
		uniqueDict = self._container._uniqueDict[self._rowIdx][self._colIdx]

		if attrs != None and len(set(attrs).difference(set(uniqueDict.keys()))):
			raise uniqueException, 'A specified attribute is not included in knoweldge base for this unique object!\n'

		if attrs == None or len(attrs) == 0 or (hasattr(attrs, "lower") and attrs.lower() == "all"):
			return uniqueDict["objects"][:]

		L = []
		for i in xrange(len(uniqueDict["objects"])):
			addThis = True
			for attr in attrs:
				if uniqueDict[attr][i] != attrs[attr]:
					addThis = False
			if addThis:
				L.append(uniqueDict["objects"][i])
		return L

	def uniqueDel(self, uniqueObj):
		# Deletes unique object uniqueObj and decrements unique count
		uniqueDict = self._container._uniqueDict[self._rowIdx][self._colIdx]
		uniqueIdx = uniqueObj._uniqueIdx

		if id(uniqueObj) != id(uniqueDict["objects"][uniqueIdx]):
			raise uniqueException, 'Unique object to delete does not match row in unique table!\n'

		for i in xrange(uniqueIdx + 1, len(uniqueDict["objects"])):
			uniqueDict["objects"][i]._uniqueIdx -= 1
		for attr in uniqueDict:
			del uniqueDict[attr][uniqueIdx]
		self._container._countsUnique[self._rowIdx, self._colIdx] -= 1			

class uniqueException(Exception):
	'''
	uniqueException
	'''


class MoleculeUniqueMeta(type):
	
	def __new__(cls, name, bases, attrs):
		attrs.update({"_container": None, "_molRowIdx": None, "_molColIdx": None})
		newClass =  super(MoleculeUniqueMeta, cls).__new__(cls, name, bases, attrs)
		_Molecule.uniqueClassRegistry[attrs["registrationId"]] = newClass
		return newClass