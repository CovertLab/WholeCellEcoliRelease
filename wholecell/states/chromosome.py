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
import random

N_BASES = 5000000
# N_CHROMOSOMES = 10

OBJ_WIDTH = 50 # width of RNA poly molecule
N_RNAP = 1000 # "anywhere from 700-2000" - nick

N_ITERS = 10 # number of testing iterations
N_CHECK = 10000 # number of regions to access
N_REMOVE = 500 # number of molecules to remove


class Molecule(object):
	# A simple object used to reference individual molecules bound to a
	# Chromosome object
	chrIndex = None
	width = OBJ_WIDTH
	boundAt = None


# TODO: handle circular chromosomes
class Chromosome(object):
	molecules = None

	empty = -1

	def __init__(self):
		# Allocate the empty chromosome object
		self.molecules = []


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
		index = self._moleculeIndex(molecule)
		molecule.boundAt = start
		self._setRange(self._range(molecule.boundAt, molecule.width), index)


	def boundMoleculeRemove(self, molecule):
		# Remove a bound molecule
		self._setRange(self._range(molecule.boundAt, molecule.width), self.empty)
		molecule.boundAt = None

		self._moleculeUnindex(molecule)


	def _moleculeIndex(self, molecule):
		if molecule.chrIndex is not None:
			return molecule.chrIndex

		else:
			try:
				newIndex = self.molecules.index(None)

			except ValueError:
				self.molecules.append(None)
				newIndex = len(self.molecules) - 1
			
			molecule.chrIndex = newIndex
			self.molecules[newIndex] = molecule

			return molecule.chrIndex


	def _moleculeUnindex(self, molecule):
		self.molecules[molecule.chrIndex] = None
		molecule.chrIndex = None


	def _setRange(self, rng, value):
		# Set some range of the chromosome representation to a value
		raise NotImplementedError()

	_range = None # Function that returns a range of values given a start and stop


class ChromosomeArray(Chromosome):
	def __init__(self):
		super(ChromosomeArray, self).__init__()
		self.chrArray = numpy.empty(N_BASES, int)
		self.chrArray[:] = self.empty


	def boundMolecules(self, start, stop):
		molInds = set(self.chrArray[start:stop]) - {self.empty}

		return {self.molecules[i] for i in molInds}


	def _setRange(self, rng, value):
		self.chrArray[rng] = value


	_range = numpy.arange


class ChromosomeDict(Chromosome):
	def __init__(self):
		super(ChromosomeDict, self).__init__()
		self.chrDict = {i:self.empty for i in self._range(N_BASES)}


	def boundMolecules(self, start, stop):
		molInds = {self.chrDict[i] for i in self._range(start, stop)} - {self.empty}

		return {self.molecules[i] for i in molInds}


	def _setRange(self, rng, value):
		for i in rng:
			self.chrDict[i] = value


	_range = xrange
		

def testClass(chromosomeClass, iters = N_ITERS):
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


def testClasses(classes = (ChromosomeArray, ChromosomeDict)):
	for cls in classes:
		timeInit, timeSetup, timeAccess, timeRemove, timeAdd, timeMove = testClass(cls)

		print ''
		print cls.__name__
		print '{:0.3f}s to initialize\n{:0.3f}s to bind {} molecules\n{:0.3f}s to check {} locations\n{:0.3f}s to remove {} molecules\n{:0.3f}s to re-bind {} molecules\n{:0.3f}s to randomly move {} molecules'.format(
			timeInit,
			timeSetup, N_RNAP,
			timeAccess, N_CHECK,
			timeRemove, N_REMOVE,
			timeAdd, N_REMOVE,
			timeMove, N_REMOVE
			)

def testSaving():
	import tables

	chromosome = ChromosomeArray()
	chromosome.setup()

	nSteps = 500

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
			expectedrows = nSteps
			)

		molecules = random.sample(chromosome.molecules, N_REMOVE)

		for t in xrange(nSteps):
			chromosome.boundMoleculeRemove(molecules[t])

			table = h5file.get_node('/', 'test')
			entry = table.row

			entry['boundMolecules'] = chromosome.chrArray+1

			entry.append()
			table.flush()

			if t % 100 == 0:
				print t


if __name__ == '__main__':
	testClasses()
