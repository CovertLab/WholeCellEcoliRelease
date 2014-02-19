'''

For the moment, this is just me playing with different imeplementations of the 
Chromosome State.

- John

'''

'''
numpy approach

Chromosomes are stored in NxM integer arrays (N bases, M multiplicity).  The 
integer is an index reference for simulation objects (specific molecules), or
-1 (nan?) if there are no objects bound.

'''

import numpy
import time

EMPTY = -1
NOT_EMPTY = 1 # until I actually do some indexing

N_BASES = 5000000
# N_CHROMOSOMES = 10

OBJ_WIDTH = 50 # width of RNA poly molecule
N_RNAP = 1000 # "anywhere from 700-2000" - nick

N_ITERS = 100 # number of testing iterations
N_CHECK = 10000 # number of regions to access

# TODO: handle circular chromosomes
class Chromosome(object):
	def __init__(self):
		# Allocate the empty chromosome object
		raise NotImplementedError()


	def setup(self):
		# Randomly attach RNAP molecules
		nBound = 0

		while nBound < N_RNAP:
			# Choose a random point and try to attach molecules

			start = numpy.random.randint(N_BASES - OBJ_WIDTH)

			if not self.boundMolecule(start, start + OBJ_WIDTH):
				self.boundMoleculeIs(NOT_EMPTY, start, OBJ_WIDTH)
				nBound += 1

	
	def boundMolecule(self, start, stop):
		# Check if a range is occupied at any point from "start" up to (not through) "stop"
		raise NotImplementedError()


	def boundMoleculeIs(self, start, width):
		# Attach a molecule to a region
		raise NotImplementedError()


class ChromosomeArray(Chromosome):
	def __init__(self):
		self.chrArray = numpy.empty(N_BASES, int)
		self.chrArray[:] = EMPTY


	def boundMolecule(self, start, stop):
		return set(self.chrArray[start:stop]) - {EMPTY}		


	def boundMoleculeIs(self, value, start, width):
		self.chrArray[start:start+width] = value


class ChromosomeDict(Chromosome):
	def __init__(self):
		self.chrDict = {i:EMPTY for i in xrange(N_BASES)}


	def boundMolecule(self, start, stop):
		return {self.chrDict[i] for i in xrange(start, stop)} - {EMPTY}


	def boundMoleculeIs(self, value, start, width):
		for i in xrange(start, start+width):
			self.chrDict[i] = value


def testClass(chromosomeClass, iters = N_ITERS):
	timeInit = 0
	timeSetup = 0
	timeAccess = 0
	
	for i in xrange(iters):
		t = time.time()
		chromosome = chromosomeClass()
		timeInit += time.time() - t

		t = time.time()
		chromosome.setup()
		timeSetup += time.time() - t

		indexes = numpy.random.randint(N_BASES - OBJ_WIDTH, size = N_CHECK)
		t = time.time()
		for ind in indexes:
			chromosome.boundMolecule(i, i + OBJ_WIDTH)
		timeAccess += time.time() - t
	
	return timeInit, timeSetup, timeAccess

for cls in [ChromosomeArray, ChromosomeDict]:
	timeInit, timeSetup, timeAccess = testClass(cls)

	print cls.__name__
	print '{:0.3f}s to initialize\n{:0.3f}s to bind {} molecules\n{:0.3f}s to check {} locations'.format(
		timeInit, timeSetup, N_RNAP, timeAccess, N_CHECK
		)

'''
ChromosomeArray
1.359s to initialize
0.928s to bind 1000 molecules
6.395s to check 10000 locations
ChromosomeDict
42.445s to initialize
1.091s to bind 1000 molecules
4.602s to check 10000 locations
'''
