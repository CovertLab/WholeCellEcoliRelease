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

import time
import random
from collections import OrderedDict

import numpy
import tables

N_BASES = 5000000
# N_CHROMOSOMES = 10

OBJ_WIDTH = 50 # width of RNA poly molecule
N_RNAP = 1000 # "anywhere from 700-2000" - nick

N_ITERS = 10 # number of testing iterations
N_CHECK = 10000 # number of regions to access
N_REMOVE = 500 # number of molecules to remove

N_STEPS = 10 # number of steps to run and save


class Molecule(object):
	# A simple object used to reference individual molecules bound to a Chromosome object
	chrIndex = None
	width = OBJ_WIDTH
	boundAt = None


# TODO: handle circular chromosomes
class Chromosome(object):
	molecules = None

	empty = -1 # TODO: concerted handling for non-molecule chromosome objects like "empty", "does-not-exist", etc.

	def __init__(self):
		self.molecules = []
		self._allocate()


	def _allocate(self):
		# Allocate the empty chromosome object
		raise NotImplementedError()


	def setup(self):
		# Randomly attach RNAP molecules
		nBound = 0

		while nBound < N_RNAP:
			# Choose a random point and try to attach molecules

			start = numpy.random.randint(N_BASES - Molecule.width)

			if not self.boundMolecules(start, start + Molecule.width):
				self.boundMoleculeIs(Molecule(), start)
				nBound += 1

	
	def boundMolecules(self, start, stop):
		# Check if a range is occupied at any point from "start" up to (not through) "stop"
		raise NotImplementedError()


	def boundMoleculeIs(self, molecule, start):
		# Attach a molecule to a region
		index = self._moleculeAssignIndex(molecule)
		molecule.boundAt = start
		self._setRange(self._range(molecule.boundAt, molecule.width), index)


	def boundMoleculeRemove(self, molecule):
		# Remove a bound molecule
		self._setRange(self._range(molecule.boundAt, molecule.width), self.empty)
		molecule.boundAt = None

		self._moleculeUnassignIndex(molecule)


	def _moleculeAssignIndex(self, molecule):
		# Assign an index to a molecule
		try:
			newIndex = self.molecules.index(None)

		except ValueError:
			self.molecules.append(None)
			newIndex = len(self.molecules) - 1
		
		molecule.chrIndex = newIndex
		self.molecules[newIndex] = molecule

		return molecule.chrIndex


	def _moleculeUnassignIndex(self, molecule):
		# Remove a molecule's index, freeing it up
		self.molecules[molecule.chrIndex] = None
		molecule.chrIndex = None


	def _setRange(self, rng, value):
		# Set some range of the chromosome representation to a value
		raise NotImplementedError()

	_range = None # Function that returns a range of values given a start and stop

	def toUnsignedArray(self):
		# Return an unsigned integer numpy array (for saving)
		# NOTE: Choosing "empty", the most common entry, to be 0, seems to
		# improve compression.
		raise NotImplementedError()


class ChromosomeArray(Chromosome):
	# numpy array implementation of the chromosome object
	def _allocate(self):
		self.chrArray = numpy.empty(N_BASES, numpy.int64)
		self.chrArray[:] = self.empty


	def boundMolecules(self, start, stop):
		molInds = set(self.chrArray[start:stop]) - {self.empty}

		return {self.molecules[i] for i in molInds}


	def _setRange(self, rng, value):
		self.chrArray[rng] = value


	def toUnsignedArray(self):
		return self.chrArray + 1


	_range = numpy.arange


class ChromosomeDict(Chromosome):
	# dictionary implementation of the chromosome object
	def _allocate(self):
		self.chrDict = {i:self.empty for i in self._range(N_BASES)}


	def boundMolecules(self, start, stop):
		molInds = {self.chrDict[i] for i in self._range(start, stop)} - {self.empty}

		return {self.molecules[i] for i in molInds}


	def _setRange(self, rng, value):
		for i in rng:
			self.chrDict[i] = value


	def toUnsignedArray(self):
		array = numpy.empty(N_BASES)

		(keys, values) = zip(*self.chrDict.items())

		array[numpy.array(keys)] = values

		array += 1

		return array


	_range = xrange


class ChromosomeOrderedDict(ChromosomeDict):
	# ordered dictionary implementation of the chromosome object
	def _allocate(self):
		self.chrDict = OrderedDict([(i, self.empty) for i in self._range(N_BASES)])


	def toUnsignedArray(self):
		return numpy.array(self.chrDict.values()) + 1
		

def testRunningClass(chromosomeClass, iters = N_ITERS):
	# Test the time needed to run a class, performing some basic operations

	timeInit = 0
	timeSetup = 0
	timeAccess = 0
	timeRemove = 0
	timeAdd = 0
	timeMove = 0
	
	for i in xrange(iters):
		t = time.time()
		chromosome = chromosomeClass()
		timeInit += time.time() - t

		t = time.time()
		chromosome.setup()
		timeSetup += time.time() - t

		indexes = numpy.random.randint(N_BASES - Molecule.width, size = N_CHECK)
		t = time.time()
		for ind in indexes:
			chromosome.boundMolecules(i, i + Molecule.width)
		timeAccess += time.time() - t

		molecules = random.sample(chromosome.molecules, N_REMOVE)
		t = time.time()
		for molecule in molecules:
			chromosome.boundMoleculeRemove(molecule)
		timeRemove += time.time() - t

		t = time.time()
		nBound = 0
		while nBound < N_REMOVE:
			start = numpy.random.randint(N_BASES - Molecule.width)

			if not chromosome.boundMolecules(start, start + Molecule.width):
				chromosome.boundMoleculeIs(Molecule(), start)
				nBound += 1	
		timeAdd += time.time() - t

		molecules = random.sample(chromosome.molecules, N_REMOVE)
		t = time.time()
		for molecule in molecules:
			newPos = min(molecule.boundAt + 40, N_BASES - molecule.width)

			if not (chromosome.boundMolecules(newPos, newPos + molecule.width) - {molecule}):
				chromosome.boundMoleculeRemove(molecule)
				chromosome.boundMoleculeIs(molecule, newPos)

		timeMove += time.time() - t
	
	return timeInit, timeSetup, timeAccess, timeRemove, timeAdd, timeMove


def testRunningClasses(classes = (ChromosomeArray, ChromosomeDict, ChromosomeOrderedDict)):
	for cls in classes:
		timeInit, timeSetup, timeAccess, timeRemove, timeAdd, timeMove = testRunningClass(cls)

		print ''
		print '{} (sum of {} runs)'.format(cls.__name__, N_ITERS)
		print '{:0.3f}s to initialize\n{:0.3f}s to bind {} molecules\n{:0.3f}s to check {} locations\n{:0.3f}s to remove {} molecules\n{:0.3f}s to re-bind {} molecules\n{:0.3f}s to randomly move {} molecules'.format(
			timeInit,
			timeSetup, N_RNAP,
			timeAccess, N_CHECK,
			timeRemove, N_REMOVE,
			timeAdd, N_REMOVE,
			timeMove, N_REMOVE,
			)


def testSaving(chromosomeClass):
	# Test the time needed to save a class while performing a simple molecule
	# unbinding operation (keeps the data from being identical between steps,
	# which could strongly bias the outcome due to compression)

	chromosome = chromosomeClass()
	chromosome.setup()

	initTime = time.time()
	with tables.open_file('test.hdf', mode = 'w', title = 'test') as h5file:
		columns = {
			'boundMolecules':tables.UInt32Col(N_BASES)
			}

		table = h5file.create_table(
			h5file.root,
			'test',
			columns,
			title = 'test',
			filters = tables.Filters(complevel = 9, complib = 'zlib'),
			expectedrows = N_STEPS
			)

		molecules = random.sample(chromosome.molecules, N_STEPS)

		for t in xrange(N_STEPS):
			chromosome.boundMoleculeRemove(molecules[t])

			table = h5file.get_node('/', 'test')
			entry = table.row

			entry['boundMolecules'] = chromosome.toUnsignedArray()

			entry.append()
			table.flush()

			# if t % 100 == 0:
			# 	print t

	return (time.time() - initTime)


def testSavingClasses(classes = (ChromosomeArray, ChromosomeDict, ChromosomeOrderedDict)):
	for cls in classes:
		timeSave = testSaving(cls)

		print ''
		print cls.__name__
		print '{:0.3f}s to save {} simulation steps'.format(timeSave, N_STEPS)


if __name__ == '__main__':
	testRunningClasses()
	testSavingClasses()

# output on my system
'''
ChromosomeArray (sum of 10 runs)
0.149s to initialize
0.318s to bind 1000 molecules
1.003s to check 10000 locations
0.014s to remove 500 molecules
0.145s to re-bind 500 molecules
0.157s to randomly move 500 molecules

ChromosomeDict (sum of 10 runs)
5.476s to initialize
0.226s to bind 1000 molecules
0.513s to check 10000 locations
0.004s to remove 500 molecules
0.105s to re-bind 500 molecules
0.102s to randomly move 500 molecules

ChromosomeOrderedDict (sum of 10 runs)
85.987s to initialize
0.266s to bind 1000 molecules
0.888s to check 10000 locations
0.005s to remove 500 molecules
0.119s to re-bind 500 molecules
0.121s to randomly move 500 molecules

ChromosomeArray
1.389s to save 10 simulation steps

ChromosomeDict
58.583s to save 10 simulation steps

ChromosomeOrderedDict
15.112s to save 10 simulation steps
'''
