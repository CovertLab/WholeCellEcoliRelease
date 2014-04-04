'''
test_chromosome_container.py

Tests for the ChromosomeContainer class.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 3/17/2014
'''

from __future__ import division

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

from wholecell.containers.chromosome_container import ChromosomeContainer

N_BASES = 1000
STRAND_MULTIPLICITY = 3

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		},
	'DNA polymerase':{
		},
	}

class Test_ChromosomeContainer(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		self.container = createContainer()


	def tearDown(self):
		pass

	# Interface tests

	# Binding and unbinding
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_bind_molecule(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 200
		forwardExtent = 5
		reverseExtent = 1

		footprint = 6
		region = [199, 200, 201, 202, 203, 204]

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '+', forwardExtent, reverseExtent)

		chromosomeIndex = mol.attr('_globalIndex') + self.container._offset
		# TODO: make the above into a private method of the container class

		# Check footprint
		self.assertEqual(
			(self.container._array[0, :] == chromosomeIndex).sum(),
			footprint
			)

		# Check location
		self.assertEqual(
			np.where(self.container._array == chromosomeIndex)[1].tolist(),
			region
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_bind_molecule_reverse(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 200
		forwardExtent = 5
		reverseExtent = 1

		footprint = 6
		region = [196, 197, 198, 199, 200, 201]

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '-', forwardExtent, reverseExtent)

		chromosomeIndex = mol.attr('_globalIndex') + self.container._offset

		# Check footprint
		self.assertEqual(
			(self.container._array[0, :] == chromosomeIndex).sum(),
			footprint
			)

		# Check location
		self.assertEqual(
			np.where(self.container._array == chromosomeIndex)[1].tolist(),
			region
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_unbind_molecule(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 200
		forwardExtent = 5
		reverseExtent = 1

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '+', forwardExtent, reverseExtent)

		self.container.moleculeLocationIsUnbound(mol)

		self.assertEqual(self.container, createContainer())


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_unbind_molecule_reverse(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 200
		forwardExtent = 5
		reverseExtent = 1

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '-', forwardExtent, reverseExtent)

		self.container.moleculeLocationIsUnbound(mol)

		self.assertEqual(self.container, createContainer())


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_bind_molecule_wrapped(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 998
		forwardExtent = 5
		reverseExtent = 1

		footprint = 6
		region = sorted([997, 998, 999, 0, 1, 2])

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '+', forwardExtent, reverseExtent)

		chromosomeIndex = mol.attr('_globalIndex') + self.container._offset

		# Check footprint
		self.assertEqual(
			(self.container._array[0, :] == chromosomeIndex).sum(),
			footprint
			)

		# Check location
		self.assertEqual(
			np.where(self.container._array == chromosomeIndex)[1].tolist(),
			region
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_bind_molecule_wrapped_reverse(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 1
		forwardExtent = 5
		reverseExtent = 1

		footprint = 6
		region = sorted([997, 998, 999, 0, 1, 2])

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '-', forwardExtent, reverseExtent)

		chromosomeIndex = mol.attr('_globalIndex') + self.container._offset

		# Check footprint
		self.assertEqual(
			(self.container._array[0, :] == chromosomeIndex).sum(),
			footprint
			)

		# Check location
		self.assertEqual(
			np.where(self.container._array == chromosomeIndex)[1].tolist(),
			region
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_unbind_molecule_wrapped(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 998
		forwardExtent = 5
		reverseExtent = 1

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '+', forwardExtent, reverseExtent)

		self.container.moleculeLocationIsUnbound(mol)

		self.assertEqual(self.container, createContainer())


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_unbind_molecule_wrapped_reverse(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 998
		forwardExtent = 5
		reverseExtent = 1

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '-', forwardExtent, reverseExtent)

		self.container.moleculeLocationIsUnbound(mol)

		self.assertEqual(self.container, createContainer())

	# Deleting
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_delete_molecule(self):
		mol = self.container.moleculeNew('DNA polymerase')

		position = 200
		forwardExtent = 5
		reverseExtent = 1

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			position, '+', forwardExtent, reverseExtent)

		self.container.moleculeDel(mol)

		self.assertEqual(self.container, createContainer())


	# Accessing bound molecules
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculesBound(self):
		positions = [100, 200, 300, 400]
		forwardExtent = 5
		reverseExtent = 1

		molecules = set()

		for position in positions:
			mol = self.container.moleculeNew('DNA polymerase')

			self.container.moleculeLocationIs(mol, self.container.rootStrand(),
				position, '+', forwardExtent, reverseExtent)

			molecules.add(mol)

		self.assertEqual(
			molecules,
			set(self.container.moleculesBound())
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculesBoundWithName(self):
		positions = [100, 200, 300, 400]
		forwardExtent = 5
		reverseExtent = 1

		molecules = set()

		for position in positions:
			mol = self.container.moleculeNew('DNA polymerase')

			self.container.moleculeLocationIs(mol, self.container.rootStrand(),
				position, '+', forwardExtent, reverseExtent)

			molecules.add(mol)

		mol = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(mol, self.container.rootStrand(),
			500, '+', forwardExtent, reverseExtent)

		self.assertNotEqual(
			molecules,
			set(self.container.moleculesBound())
			)

		self.assertEqual(
			molecules,
			set(self.container.moleculesBoundWithName('DNA polymerase'))
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculeBoundAtPosition(self):
		positions = [100, 200, 300, 400]
		forwardExtent = 5
		reverseExtent = 1

		molecules = []

		for position in positions:
			mol = self.container.moleculeNew('DNA polymerase')

			self.container.moleculeLocationIs(mol, self.container.rootStrand(),
				position, '+', forwardExtent, reverseExtent)

			molecules.append(mol)

		self.assertEqual(
			molecules[1],
			self.container.moleculeBoundAtPosition(
				self.container.rootStrand(), 200
				)
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculeBoundOnFork(self):
		startPosition = 100
		stopPosition = 110

		forkStart, forkStop = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		forwardExtent = 5
		reverseExtent = 1

		mol = self.container.moleculeNew('DNA polymerase')

		self.container.moleculeLocationIsFork(mol, forkStart, forwardExtent,
			reverseExtent)

		self.assertEqual(
			mol,
			self.container.moleculeBoundOnFork(forkStart)
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculesBoundOverExtent(self):
		positions = [100, 200, 300, 400]
		forwardExtent = 5
		reverseExtent = 1

		molecules = []

		for position in positions:
			mol = self.container.moleculeNew('DNA polymerase')

			self.container.moleculeLocationIs(mol, self.container.rootStrand(),
				position, '+', forwardExtent, reverseExtent)

			molecules.append(mol)

		self.assertEqual(
			set(molecules[:2]),
			set(self.container.moleculesBoundOverExtent(
				self.container.rootStrand(), 0, '+',
				200, 0
				))
			)


	# Forks
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_new_forks(self):
		startPosition = 100
		stopPosition = 110

		forks = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		footprint = 11
		region = [100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110]

		self.assertEqual(
			(self.container._array[0, :] == self.container._inactive).sum(),
			footprint - 2 # 2 spots occupied by forks
			)

		self.assertEqual(
			np.where(self.container._array[0, :] != self.container._empty)[0].tolist(),
			region
			)

		self.assertEqual(
			(self.container._array[1, :] == self.container._empty).sum(),
			footprint
			)

		self.assertEqual(
			(self.container._array[2, :] == self.container._empty).sum(),
			footprint
			)

		self.assertEqual(
			forks[0].attr('_chromPosition'),
			startPosition
			)

		self.assertEqual(
			forks[1].attr('_chromPosition'),
			stopPosition
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_new_forks_wrapped(self):
		startPosition = 995
		stopPosition = 5

		forks = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		footprint = 11
		region = sorted([995, 996, 997, 998, 999, 0, 1, 2, 3, 4, 5])

		self.assertEqual(
			(self.container._array[0, :] == self.container._inactive).sum(),
			footprint - 2 # 2 spots occupied by forks
			)

		self.assertEqual(
			np.where(self.container._array[0, :] != self.container._empty)[0].tolist(),
			region
			)

		self.assertEqual(
			(self.container._array[1, :] == self.container._empty).sum(),
			footprint
			)

		self.assertEqual(
			(self.container._array[2, :] == self.container._empty).sum(),
			footprint
			)

		self.assertEqual(
			forks[0].attr('_chromPosition'),
			startPosition
			)

		self.assertEqual(
			forks[1].attr('_chromPosition'),
			stopPosition
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_extend_forks(self):
		startPosition = 100
		stopPosition = 110

		[startFork, stopFork] = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		newStart = self.container.forkExtend(startFork, 2)
		newStop = self.container.forkExtend(stopFork, 2)

		footprint = 15
		region = [98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
			110, 111, 112]

		self.assertEqual(
			(self.container._array[0, :] == self.container._inactive).sum(),
			footprint - 2 # 2 spots occupied by forks
			)

		self.assertEqual(
			np.where(self.container._array[0, :] != self.container._empty)[0].tolist(),
			region
			)

		self.assertEqual(newStart, 98)
		self.assertEqual(newStop, 112)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_extend_forks_wrapped(self):
		startPosition = 995
		stopPosition = 5

		[startFork, stopFork] = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		newStart = self.container.forkExtend(startFork, 2)
		newStop = self.container.forkExtend(stopFork, 2)

		footprint = 15
		region = sorted([993, 994, 995, 996, 997, 998, 999, 0, 1, 2, 3, 4, 5, 6, 7])

		self.assertEqual(
			(self.container._array[0, :] == self.container._inactive).sum(),
			footprint - 2 # 2 spots occupied by forks
			)

		self.assertEqual(
			np.where(self.container._array[0, :] != self.container._empty)[0].tolist(),
			region
			)

		self.assertEqual(newStart, 993)
		self.assertEqual(newStop, 7)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_extend_fork_over_big_edge(self):
		startPosition = 994
		stopPosition = 995

		[startFork, stopFork] = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		newStop = self.container.forkExtend(stopFork, 10)

		footprint = 12
		region = sorted([994, 995, 996, 997, 998, 999, 0, 1, 2, 3, 4, 5])

		self.assertEqual(
			(self.container._array[0, :] == self.container._inactive).sum(),
			footprint - 2 # 2 spots occupied by forks
			)

		self.assertEqual(
			np.where(self.container._array[0, :] != self.container._empty)[0].tolist(),
			region
			)

		self.assertEqual(newStop, 5)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_extend_fork_over_small_edge(self):
		startPosition = 1
		stopPosition = 3

		[startFork, stopFork] = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		newStart = self.container.forkExtend(startFork, 10)

		footprint = 13
		region = sorted([991, 992, 993, 994, 995, 996, 997, 998, 999, 0, 1, 2, 3])

		self.assertEqual(
			(self.container._array[0, :] == self.container._inactive).sum(),
			footprint - 2 # 2 spots occupied by forks
			)

		self.assertEqual(
			np.where(self.container._array[0, :] != self.container._empty)[0].tolist(),
			region
			)

		self.assertEqual(newStart, 991)


	# Molecules on forks
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_bind_molecule_to_fork(self):
		mol = self.container.moleculeNew('DNA polymerase')

		startPosition = 100
		stopPosition = 110

		[forkStart, forkStop] = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		forwardExtent = 5
		reverseExtent = 1

		footprintParent = 5
		regionParent = [111, 112, 113, 114, 115]

		footprintChild = 1
		regionChild = [110]

		self.container.moleculeLocationIsFork(mol, forkStop,
			forwardExtent, reverseExtent)

		chromosomeIndex = mol.attr('_globalIndex') + self.container._offset

		# Check footprint
		self.assertEqual(
			(self.container._array[0, :] == chromosomeIndex).sum(),
			footprintParent
			)

		self.assertEqual(
			(self.container._array[1, :] == chromosomeIndex).sum(),
			footprintChild
			)

		self.assertEqual(
			(self.container._array[2, :] == chromosomeIndex).sum(),
			footprintChild
			)

		# Check location
		self.assertEqual(
			np.where(self.container._array[0, :] == chromosomeIndex)[0].tolist(),
			regionParent
			)

		self.assertEqual(
			np.where(self.container._array[1, :] == chromosomeIndex)[0].tolist(),
			regionChild
			)

		self.assertEqual(
			np.where(self.container._array[2, :] == chromosomeIndex)[0].tolist(),
			regionChild
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_bind_molecule_to_fork_reverse(self):
		mol = self.container.moleculeNew('DNA polymerase')

		startPosition = 100
		stopPosition = 110

		[forkStart, forkStop] = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		forwardExtent = 5
		reverseExtent = 1

		footprintParent = 5
		regionParent = [95, 96, 97, 98, 99]

		footprintChild = 1
		regionChild = [100]

		self.container.moleculeLocationIsFork(mol, forkStart,
			forwardExtent, reverseExtent)

		chromosomeIndex = mol.attr('_globalIndex') + self.container._offset

		# Check footprint
		self.assertEqual(
			(self.container._array[0, :] == chromosomeIndex).sum(),
			footprintParent
			)

		self.assertEqual(
			(self.container._array[1, :] == chromosomeIndex).sum(),
			footprintChild
			)

		self.assertEqual(
			(self.container._array[2, :] == chromosomeIndex).sum(),
			footprintChild
			)

		# Check location
		self.assertEqual(
			np.where(self.container._array[0, :] == chromosomeIndex)[0].tolist(),
			regionParent
			)

		self.assertEqual(
			np.where(self.container._array[1, :] == chromosomeIndex)[0].tolist(),
			regionChild
			)

		self.assertEqual(
			np.where(self.container._array[2, :] == chromosomeIndex)[0].tolist(),
			regionChild
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_unbind_molecule_on_fork(self):
		startPosition = 100
		stopPosition = 110

		[forkStart, forkStop] = self.container.divideRegion(
			self.container.rootStrand(),
			startPosition, stopPosition
			)

		mol = self.container.moleculeNew('DNA polymerase')

		forwardExtent = 5
		reverseExtent = 1

		footprintParent = 4 # molecule origin is occupied by fork on parent
		regionParent = [111, 112, 113, 114]

		footprintChild = 2
		regionChild = [109, 110]

		self.container.moleculeLocationIsFork(mol, forkStop,
			forwardExtent, reverseExtent)

		self.container.moleculeLocationIsUnbound(mol)

		newContainer = createContainer()

		newContainer.divideRegion(
			newContainer.rootStrand(),
			startPosition, stopPosition
			)

		self.assertEqual(self.container, newContainer)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculeFootprint(self):
		molecule = self.container.moleculeNew('RNA polymerase')

		strand = self.container.rootStrand()
		position = 150
		direction = '+'
		forwardExtent = 20
		reverseExtent = 10

		self.container.moleculeLocationIs(molecule, strand, position, 
			direction, forwardExtent, reverseExtent)

		footprint = self.container.moleculeFootprint(molecule)

		self.assertEqual(
			set(footprint),
			set(np.arange(140, 170))
			)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculesBoundPastFork(self):
		strand = self.container.rootStrand()

		fork = self.container.divideRegion(strand, 100, 200)[1]

		positions = [210, 230, 280]
		direction = '+'
		forwardExtent = 5
		reverseExtent = 5

		molecules = []

		for position in positions:
			molecule = self.container.moleculeNew('RNA polymerase')

			self.container.moleculeLocationIs(molecule, strand, position,
				direction, forwardExtent, reverseExtent)

			molecules.append(molecule)

		self.assertEqual(
			set(self.container.moleculesBoundPastFork(fork, 4)),
			set()
			)

		self.assertEqual(
			set(self.container.moleculesBoundPastFork(fork, 10)),
			set(molecules[:1])
			)

		self.assertEqual(
			set(self.container.moleculesBoundPastFork(fork, 90)),
			set(molecules[:3])
			)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_regionsNearForks_simple(self):
		strand = self.container.rootStrand()

		self.container.divideRegion(strand, 100, 200)

		forwardExtent = 4
		reverseExtent = 2

		regionsParent, regionsChildA, regionsChildB = self.container.regionsNearForks(
			forwardExtent, reverseExtent, False)

		indexesFork1 = {96, 97, 98, 99, 100, 101}
		indexesFork2 = {199, 200, 201, 202, 203, 204}

		for region in regionsParent:
			self.assertEqual(region.strand(), 0)

			indexes = region.indexes()

			if indexes[0] < 150:
				self.assertEqual(
					set(indexes),
					indexesFork1
					)

			else:
				self.assertEqual(
					set(indexes),
					indexesFork2
					)

		for region in regionsChildA:
			self.assertEqual(region.strand(), 1)

			indexes = region.indexes()

			if indexes[0] < 150:
				self.assertEqual(
					set(indexes),
					indexesFork1
					)

			else:
				self.assertEqual(
					set(indexes),
					indexesFork2
					)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_regionsNearForks_minimal_extents(self):
		strand = self.container.rootStrand()

		forks = self.container.divideRegion(strand, 100, 200)

		molExtentForward = 4
		molExtentReverse = 2

		for fork in forks:
			molecule = self.container.moleculeNew('DNA polymerase')
			self.container.moleculeLocationIsFork(molecule, fork,
				molExtentForward, molExtentReverse)

		regionsParent, regionsChildA, regionsChildB = self.container.regionsNearForks(
			1, 1, False)

		indexesFork1 = {96, 97, 98, 99, 100, 101}
		indexesFork2 = {199, 200, 201, 202, 203, 204}

		for region in regionsParent:
			self.assertEqual(region.strand(), 0)

			indexes = region.indexes()

			if indexes[0] < 150:
				self.assertEqual(
					set(indexes),
					indexesFork1
					)

			else:
				self.assertEqual(
					set(indexes),
					indexesFork2
					)

		for region in regionsChildA:
			self.assertEqual(region.strand(), 1)

			indexes = region.indexes()

			if indexes[0] < 150:
				self.assertEqual(
					set(indexes),
					indexesFork1
					)

			else:
				self.assertEqual(
					set(indexes),
					indexesFork2
					)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_regionsNearForks_include_ends(self):
		strand = self.container.rootStrand()

		self.container.divideRegion(strand, 100, 200)

		molecule1 = self.container.moleculeNew('RNA polymerase')
		molecule2 = self.container.moleculeNew('RNA polymerase')
		molecule3 = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			molecule1,
			strand, 90, '-',
			10, 6
			)

		self.container.moleculeLocationIs(
			molecule2,
			strand + 'A', 190, '+',
			10, 6
			)

		self.container.moleculeLocationIs(
			molecule3,
			strand, 210, '+',
			10, 6
			)

		forwardExtent = 4
		reverseExtent = 2

		regionsParent, regionsChildA, regionsChildB = self.container.regionsNearForks(
			forwardExtent, reverseExtent, True)

		indexesFork1 = set(range(81, 102))
		indexesFork2 = set(range(199, 220))

		for region in regionsParent:
			self.assertEqual(region.strand(), 0)

			indexes = region.indexes()

			if indexes[0] < 150:
				self.assertEqual(
					set(indexes),
					indexesFork1
					)

			else:
				self.assertEqual(
					set(indexes),
					indexesFork2
					)

		for region in regionsChildA:
			self.assertEqual(region.strand(), 1)

			indexes = region.indexes()

			if indexes[0] < 150:
				self.assertEqual(
					set(indexes),
					indexesFork1
					)

			else:
				self.assertEqual(
					set(indexes),
					indexesFork2
					)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_regionsNearMolecules_simple(self):
		strand = self.container.rootStrand()

		molecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			molecule,
			strand, 100, '+',
			4, 2
			)

		forwardExtent = 10
		reverseExtent = 5

		regions = self.container.regionsNearMolecules([molecule],
			forwardExtent, reverseExtent, False)

		self.assertEqual(len(regions), 1)

		(region,) = regions

		self.assertEqual(
			set(region.indexes()),
			set(range(95, 110))
			)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_regionsNearMolecules_minimal_extents(self):
		strand = self.container.rootStrand()

		molecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			molecule,
			strand, 100, '+',
			4, 2
			)

		forwardExtent = 0
		reverseExtent = 0

		regions = self.container.regionsNearMolecules([molecule],
			forwardExtent, reverseExtent, False)

		self.assertEqual(len(regions), 1)

		(region,) = regions

		self.assertEqual(
			set(region.indexes()),
			set(range(98, 104))
			)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_regionsNearMolecules_include_ends(self):
		strand = self.container.rootStrand()

		molecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			molecule,
			strand, 100, '+',
			4, 2
			)

		endMolecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			endMolecule,
			strand, 110, '+',
			4, 2
			)

		forwardExtent = 10
		reverseExtent = 5

		regions = self.container.regionsNearMolecules([molecule],
			forwardExtent, reverseExtent, True)

		self.assertEqual(len(regions), 1)

		(region,) = regions

		self.assertEqual(
			set(region.indexes()),
			set(range(95, 114))
			)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculeInRegionSet(self):
		strand = self.container.rootStrand()

		molecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			molecule,
			strand, 100, '+',
			4, 2
			)

		endMolecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			endMolecule,
			strand, 110, '+',
			4, 2
			)

		otherMolecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			otherMolecule,
			strand, 200, '+',
			4, 2
			)

		forwardExtent = 10
		reverseExtent = 5

		regions = self.container.regionsNearMolecules([molecule],
			forwardExtent, reverseExtent, True)

		self.assertTrue(
			self.container.moleculeInRegionSet(molecule, regions)
			)

		self.assertTrue(
			self.container.moleculeInRegionSet(endMolecule, regions)
			)

		self.assertFalse(
			self.container.moleculeInRegionSet(otherMolecule, regions)
			)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculesInRegion(self):
		strand = self.container.rootStrand()

		molecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			molecule,
			strand, 100, '+',
			4, 2
			)

		endMolecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			endMolecule,
			strand, 110, '+',
			4, 2
			)

		otherMolecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			otherMolecule,
			strand, 200, '+',
			4, 2
			)

		forwardExtent = 10
		reverseExtent = 5

		regions = self.container.regionsNearMolecules([molecule],
			forwardExtent, reverseExtent, True)

		(region,) = regions

		self.assertEqual(
			set(self.container.moleculesInRegion(region)),
			{molecule, endMolecule}
			)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculesInRegionSet(self):
		strand = self.container.rootStrand()

		positions = [20, 40, 60, 80, 100]
		molecules = []

		for position in positions:
			molecule = self.container.moleculeNew('RNA polymerase')

			self.container.moleculeLocationIs(
				molecule,
				strand, position, '+',
				4, 2
				)

			molecules.append(molecule)

		otherMolecule = self.container.moleculeNew('RNA polymerase')

		self.container.moleculeLocationIs(
			otherMolecule,
			strand, 120, '+',
			4, 2
			)

		forwardExtent = 10
		reverseExtent = 5

		regions = self.container.regionsNearMolecules(molecules,
			forwardExtent, reverseExtent, False)

		moleculesInRegionSet = self.container.moleculesInRegionSet(regions)

		self.assertEqual(
			set(moleculesInRegionSet),
			set(molecules)
			)


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_forksInRegion(self):
		pass


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_forkInRegion(self):
		pass


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_forkInRegionSet(self):
		pass


	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_maximumExtentPastFork(self):
		pass

	# TODO: specific _ChromosomeRegion[Set] tests


def createContainer():
	container = ChromosomeContainer(N_BASES, STRAND_MULTIPLICITY,
		MOLECULE_ATTRIBUTES)

	return container
