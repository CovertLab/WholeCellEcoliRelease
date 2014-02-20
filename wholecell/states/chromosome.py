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

EMPTY = -1

N_BASES = 5000000
# N_CHROMOSOMES = 10

OBJ_WIDTH = 50 # width of RNA poly molecule
N_RNAP = 1000 # "anywhere from 700-2000" - nick

N_ITERS = 100 # number of testing iterations
N_CHECK = 10000 # number of regions to access
N_REMOVE = 500 # number of molecules to remove


class Molecule(object):
	# A simple object used to reference individual molecules bound to a
	# Chromosome object
	chrIndex = None
	boundRange = None


# TODO: handle circular chromosomes
# TODO: push reindexing to a cleanup routine
class Chromosome(object):
	molecules = None
	def __init__(self):
		# Allocate the empty chromosome object
		self.molecules = []


	def setup(self):
		# Randomly attach RNAP molecules
		nBound = 0

		while nBound < N_RNAP:
			# Choose a random point and try to attach molecules

			start = numpy.random.randint(N_BASES - OBJ_WIDTH)

			if not self.boundMolecules(start, start + OBJ_WIDTH):
				self.boundMoleculeIs(Molecule(), start, OBJ_WIDTH)
				nBound += 1

	
	def boundMolecules(self, start, stop):
		# Check if a range is occupied at any point from "start" up to (not through) "stop"
		raise NotImplementedError()


	def boundMoleculeIs(self, molecule, start, width):
		# Attach a molecule to a region
		raise NotImplementedError()


	def boundMoleculeRemove(self, molecule):
		# Remove a molecule from the chromosome
		raise NotImplementedError()


	# def boundMoleculeMove(self, molecule, )


	def _moleculeIndex(self, molecule):
		if molecule.chrIndex is not None:
			return molecule.chrIndex

		else:
			self.molecules.append(molecule)
			return len(self.molecules) - 1


class ChromosomeArray(Chromosome):
	def __init__(self):
		super(ChromosomeArray, self).__init__()
		self.chrArray = numpy.empty(N_BASES, int)
		self.chrArray[:] = EMPTY


	def boundMolecules(self, start, stop):
		molInds = set(self.chrArray[start:stop]) - {EMPTY}

		return {self.molecules[i] for i in molInds}


	def boundMoleculeIs(self, molecule, start, width):
		index = self._moleculeIndex(molecule)
		self.chrArray[start:start+width] = index

		molecule.chrIndex = index
		molecule.boundRange = numpy.arange(start, start+width)


	def boundMoleculeRemove(self, molecule):
		self.chrArray[molecule.boundRange] = EMPTY

		for otherMol in self.molecules[molecule.chrIndex+1:]:
			otherMol.chrIndex -= 1
			self.chrArray[otherMol.boundRange] = otherMol.chrIndex

		del self.molecules[molecule.chrIndex]
		molecule.chrIndex = None
		molecule.boundRange = None


class ChromosomeDict(Chromosome):
	def __init__(self):
		super(ChromosomeDict, self).__init__()
		self.chrDict = {i:EMPTY for i in xrange(N_BASES)}


	def boundMolecules(self, start, stop):
		molInds = {self.chrDict[i] for i in xrange(start, stop)} - {EMPTY}

		return {self.molecules[i] for i in molInds}


	def boundMoleculeIs(self, molecule, start, width):
		index = self._moleculeIndex(molecule)
		for i in xrange(start, start+width):
			self.chrDict[i] = index

		molecule.chrIndex = index
		molecule.boundRange = xrange(start, start+width)


	def boundMoleculeRemove(self, molecule):
		for i in molecule.boundRange:
			self.chrDict[i] = EMPTY

		for otherMol in self.molecules[molecule.chrIndex+1:]:
			otherMol.chrIndex -= 1
			for i in molecule.boundRange:
				self.chrDict[i] = otherMol.chrIndex

		del self.molecules[molecule.chrIndex]
		molecule.chrIndex = None
		molecule.boundRange = None
		

def testClass(chromosomeClass, iters = N_ITERS):
	timeInit = 0
	timeSetup = 0
	timeAccess = 0
	timeRemove = 0
	
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
			chromosome.boundMolecules(i, i + OBJ_WIDTH)
		timeAccess += time.time() - t

		molecules = random.sample(chromosome.molecules, N_REMOVE)
		t = time.time()
		for molecule in molecules:
			chromosome.boundMoleculeRemove(molecule)
		timeRemove += time.time() - t
	
	return timeInit, timeSetup, timeAccess, timeRemove


def testClasses():
	for cls in [ChromosomeArray, ChromosomeDict]:
		timeInit, timeSetup, timeAccess, timeRemove = testClass(cls)

		print ''
		print cls.__name__
		print '{:0.3f}s to initialize\n{:0.3f}s to bind {} molecules\n{:0.3f}s to check {} locations\n{:0.3f}s to remove {} molecules'.format(
			timeInit, timeSetup, N_RNAP, timeAccess, N_CHECK, timeRemove, N_REMOVE
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


