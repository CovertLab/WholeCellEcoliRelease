""" MoleculesContainer """

import numpy
import types
import sys
import time

# TODO: Change 'object' to wholecell.sim.state.State
class MoleculesContainer(object):

	def __init__(self):
		self._countsBulk = None
		self._countsUnique = None
		self._massSingle = None
		self._dMass = None
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
			self._molecules[(wid, comp)] = Molecule(self, self._widIdx[wid], self._cIdx[comp], wid)
		return self._molecules[(wid, comp)]

class Molecule(object):
	uniqueClassRegistry = {}
	def __init__(self, container, rowIdx, colIdx, wid):
		self._container = container
		self._rowIdx = rowIdx
		self._colIdx = colIdx
		self._wid = wid

		if len(self._container._uniqueDict[self._rowIdx][self._colIdx]) == 0:
			return

		def uniqueInit(self, uniqueIdx):
			self._uniqueIdx = uniqueIdx

		if self._wid in self.uniqueClassRegistry:
			self.MoleculeUnique = self.uniqueClassRegistry[self._wid]
			self.MoleculeUnique._container = self._container
			self.MoleculeUnique._molRowIdx = self._rowIdx
			self.MoleculeUnique._molColIdx = self._colIdx
			if type(self.MoleculeUnique.__init__) != type(self.__init__):
				setattr(self.MoleculeUnique, "__init__", uniqueInit)
			return

		uniqueClassDefDict = {}
		uniqueClassDefDict["_container"] = self._container
		uniqueClassDefDict["_molRowIdx"] = self._rowIdx
		uniqueClassDefDict["_molColIdx"] = self._colIdx



		uniqueClassDefDict["__init__"] = uniqueInit

		def makeGetter(attr):
			def attrGetter(self):
				molCont, uniqueIdx = self._container, self._uniqueIdx
				molRowIdx, molColIdx = self._molRowIdx, self._molColIdx

				return molCont._uniqueDict[molRowIdx][molColIdx][attr][uniqueIdx]
			return attrGetter
		
		def makeSetter(attr):
			def attrSetter(self, newVal):
				molCont, uniqueIdx = self._container, self._uniqueIdx
				molRowIdx, molColIdx = self._molRowIdx, self._molColIdx

				molCont._uniqueDict[molRowIdx][molColIdx][attr][uniqueIdx] = newVal
			return attrSetter

		for attr in self._container._uniqueDict[self._rowIdx][self._colIdx]:

			uniqueClassDefDict[attr] = makeGetter(attr)
			uniqueClassDefDict[attr + "Is"] = makeSetter(attr)

		self.MoleculeUnique = type("MoleculeUnique", (), uniqueClassDefDict)

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
			raise uniqueException, 'Attribute not included in knoweldge base for this unique object!\n'

		for attr in uniqueDict:
			if attrs != None and attr in attrs:
				uniqueDict[attr].append(attrs[attr])
			else:
				uniqueDict[attr].append(None)
		uniqueIdx = len(uniqueDict["objects"]) - 1
		uniqueDict["objects"][uniqueIdx] = self.MoleculeUnique(uniqueIdx)
		self._container._countsUnique[self._rowIdx, self._colIdx] += 1
		return uniqueDict["objects"][uniqueIdx]

	def uniquesWithAttrs(self, attrs = None):
		uniqueDict = self._container._uniqueDict[self._rowIdx][self._colIdx]

		if attrs == None or len(attrs) == 0 or (hasattr(attrs, "lower") and attrs.lower() == "all"):
			return uniqueDict["objects"][:]

		L = []
		for i in xrange(len(uniqueDict["objects"])):
			for attr in attrs:
				if uniqueDict[attr][i] == attrs[attr]:
					L.append(uniqueDict["objects"][i])
		return L

	def uniqueDel(self, uniqueObj):
		uniqueDict = self._container._uniqueDict[self._rowIdx][self._colIdx]
		uniqueIdx = uniqueObj._uniqueIdx
		for i in xrange(uniqueIdx + 1, len(uniqueDict["objects"])):
			uniqueDict["objects"][i]._uniqueIdx -= 1
		for attr in uniqueDict:
			del uniqueDict[attr][uniqueIdx]
		self._container._countsUnique[self._rowIdx, self._colIdx] -= 1			

class MoleculeUniqueMeta(type):
	
	def __new__(cls, name, bases, attrs):
		attrs.update({"_container": None, "_molRowIdx": None, "_molColIdx": None})
		newClass =  super(MoleculeUniqueMeta, cls).__new__(cls, name, bases, attrs)
		Molecule.uniqueClassRegistry[attrs["registrationId"]] = newClass
		return newClass

class uniqueException(Exception):
	'''
	uniqueException
	'''

###
# class enz4Unique(object):
# 	registrationId = "enz4"
# 	__metaclass__ = MoleculeUniqueMeta

# #	def __init__(self, uniqueIdx):
# #		self._uniqueIdx = uniqueIdx

# 	def attr1(self):
# 		return self._container._uniqueDict[self._molRowIdx][self._molColIdx]["attr4_1"][self._uniqueIdx]

# class rnas(object):
# 	registrationType = "rnas"
# 	__metaclass__ = MoleculeUniqueMeta

##
# class Track(object):

# 	def __init__(self):
# 		self._d = {}
# 		self._maxKey = 0
		
# 	def __getitem__(self, key):
# 		if type(key) == int:
# 			if key < 0:
# 				raise IndexError, "Key must be non-negative."
# 			return self._d.get(key, None)
# 		elif type(key) == slice or type(key) == list or type(key) == tuple or type(key) == xrange or type(key) == numpy.ndarray:
# 			if type(key) == slice:
# 				seq = xrange(*key.indices(self._maxKey + 1))
# 			else:
# 				seq = key

# 			if seq[0] < 0 or seq[-1] < 0:
# 				raise IndexError, "Keys must be non-negative."

# 			d = {}
# 			for i in seq:
# 				if i in self._d:
# 					d[i] = self._d[i]
# 			return d
# 		else:
# 			raise TypeError, "Invalid argument type."

# 	def __setitem__(self, key, value):
# 		if type(key) == int:
# 			if key < 0:
# 				raise IndexError, "Key must be non-negative."
# 			if value == None:
# 				del self[key]
# 			else:
# 				self._maxKey = max(self._maxKey, key)
# 				self._d[key] = value
# 		elif type(key) == slice or type(key) == list or type(key) == tuple or type(key) == xrange or type(key) == numpy.ndarray:
# 			if type(key) == slice:
# 				seq = xrange(*key.indices(sys.maxint))
# 			else:
# 				seq = key

# 			if type(value) != list and type(value) != tuple and type(value) != xrange and type(value) != numpy.ndarray:
# 				valueIter = (value for x in xrange(len(seq)))
# 			else:
# 				if len(seq) != len(value):
# 					raise ValueError, "Shape mismatch."
# 				valueIter = iter(value)

# 			for iSeq in seq:
# 				self[iSeq] = valueIter.next()
# 				print self[iSeq]
# 		else:
# 			raise TypeError, "Invalid argument type."

# 	def __delitem__(self, key):
# 		if type(key) == int:
# 			if key < 0:
# 				raise IndexError, "Key must be non-negative."
# 			del self._d[key]
# 			if key == self._maxKey:
# 				self._maxKey = reduce(lambda x, y: max(x, y), self._d.iterkeys())
# 		elif type(key) == slice or type(key) == list or type(key) == tuple or type(key) == xrange or type(key) == numpy.ndarray:
# 			if type(key) == slice:
# 				seq = xrange(*key.indices(self._maxKey))
# 			else:
# 				seq = key

# 			d = {}
# 			for i in seq:
# 				del self[i]
# 			return d
# 		else:
# 			raise TypeError, "Invalid argument type."

# 	def __repr__(self):
# 		return repr(self._d)

# class RNA(object):

# 	def __init__(self, rnaId):
# 		self._rnaId = rnaId
# 		# TODO: Add uniqueIdx
# 		self.tracks = []
# 		self.tracksIdxs = {
# 			"sequence": [0],
# 			"bound": [1],
# 			"secStruct": [2],
# 			"forPolypeptide": [3],
# 		}
# 		[self.tracks.append(Track()) for key in self.tracksIdxs for elem in self.tracksIdxs[key]]

# class RNAP(object):

# 	ftpt = 50
# 	velocity = 50

# 	def __init__(self):
# 		# TODO: Add uniqueIdx
# 		self._state = "F"
# 		self._rna = None
# 		self._chromLocs = -1 * numpy.ones(self.ftpt, dtype = int)
# 		self._rnaLocs = [-1]

# 	def state(self): return self._state
# 	def stateIs(self, val): self._state = val
# 	def rna(self): return self._rna
# 	def rnaIs(self, val): self._rna = val
# 	def chromLocs(self, i): return self._chromLocs[i]
# 	def chromLocsLen(self): return len(self._chromLocs)
# 	def chromLocsIs(self, i, val): self._chromLocs[i] = val
# 	def chromLocsIncr(self, i, val = None):
# 		if val == None:
# 			self._chromLocs = (self._chromLocs + i) % 5e6
# 		else:
# 			self._chromLocs[i] = (self._chromLocs[i] + val) % 5e6
# 	def rnaLocs(self, i): return self._rnaLocs[i]
# 	def rnaLocsLen(self): return len(self._rnaLocs)
# 	def rnaLocsIs(self, i, val): self._rnaLocs[i] = val
# 	def rnaLocsIncr(self, i, val): self._rnaLocs[i] += val

# numpy.random.seed(1)
# chrLen = 5e6
# dna = "".join(["ACGU"[x] for x in numpy.random.choice(numpy.arange(4), size = chrLen)])
# rnaps = [RNAP() for x in xrange(500)]
# rnas = []

# before = time.time()
# for t in xrange(10):
# 	freeRnaps = [x for x in rnaps if x.state() == "F"]
# 	for rnap in freeRnaps:
# 		if numpy.random.uniform() > 0.5: continue
# 		pos = numpy.random.random_integers(0, chrLen)
# 		[rnap.chromLocsIs(x, (x + pos) % chrLen) for x in xrange(rnap.chromLocsLen())]
# 		rnap.stateIs("B")
# 		rnap.rnaIs(RNA(str(len(rnas) + 1)))
# 		rnap.rnaLocsIs(0, 0)
# 		rnas.append(rnap.rna())

# 	ntCounts = numpy.random.poisson(1e4, 4)
# 	ntLookup = {"A": 0, "C": 1, "G": 2, "U": 3}
# 	boundRnaps = [x for x in rnaps if x.state() == "B"]
# 	rnapsDone = numpy.zeros(len(boundRnaps))
# 	rnapsStuck = numpy.zeros(len(boundRnaps))
# 	# while not numpy.all(rnapsDone) and not numpy.all(rnapsStuck):
# 	for rnapIdx, rnap in enumerate(boundRnaps):
# 		if rnap.state() == "F": continue
# 		for relPos in xrange(rnap.velocity):
# 			meanPos = int(numpy.mean((rnap.chromLocs(0), rnap.chromLocs(-1))))
# 			rna = rnap.rna()
# 			for x in rna.tracksIdxs["sequence"]: rna.tracks[x][rnap.rnaLocs(0)] = dna[meanPos]
# 			# [rnap.chromLocsIncr(x, 1) for x in xrange(rnap.chromLocsLen())]
# 			rnap.chromLocsIncr(1)
# 			rnap.rnaLocsIncr(0, 1)
# 			if numpy.random.uniform() < 0.001:
# 				rnap.rnaIs(None)
# 				[rnap.chromLocsIs(x, -1) for x in xrange(rnap.chromLocsLen())]
# 				rnap.stateIs("F")
# 				break

# after = time.time()
# print "Elapsed: %0.3f seconds" % (after - before)
# class TrackGroup(object):

# 	def __init__(self):
# 		self._d = {}
# 		self._maxKey = 0

# 	def __getitem__(self, key):
# 		if type(key) == int:
# 			if key < 0:
# 				raise IndexError, "Key must be non-negative."
# 			return self._d[key]
# 		elif type(key) == slice:

# 			if key.start == None and key.step == None and key.stop == None:
# 				return [x for x in self._d.itervalues()]
# 		elif type(key) == tuple and len(key) == 2:
# 			pass
# 		else:
# 			raise TypeError, "Invalid argument type."

# 	def __setitem__(self, key, value):
# 		if type(key) == int:
# 			if key < 0:
# 				raise IndexError, "Key must be non-negative."
# 			if value == None:
# 				del self[key]
# 			else:
# 				self._maxKey = max(self._maxKey, key)
# 				self._d[key] = value
# 		else:
# 			raise TypeError, "Invalid argument type."

# 	def __delitem__(self, key):
# 		if type(key) == int:
# 			if key < 0:
# 				raise IndexError, "Key must be non-negative."
# 			del self._d[key]
# 			if key == self._maxKey:
# 				self._maxKey = reduce(lambda x, y: max(x, y), self._d.iterkeys())
# 		else:
# 			raise TypeError, "Invalid argument type."


# kb = type("", (), {})()
# kb.molecules = [{"id": "enz1", "mass": 1.0, "uniqueAttrs": None}, {"id": "enz2", "mass": 2.0, "uniqueAttrs": None}, {"id": "enz3", "mass": 3.0, "uniqueAttrs": ["attr1", "attr2", "attr3"]}, {"id": "enz4", "mass": 4.0, "uniqueAttrs": ["attr4_1", "attr4_2"]}, {"id": "enz5", "mass": 5.0, "uniqueAttrs": None}]
# kb.compartments = [{"id": "c"}, {"id": "e"}, {"id": "m"}]

# mc = MoleculesContainer()
# mc.initialize(kb)
# mol = mc.molecule("enz3", "c")
