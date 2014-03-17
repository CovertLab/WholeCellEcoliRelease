'''
test_chromosome_bound_molecules_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 3/17/2014
'''

from __future__ import division

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

from wholecell.utils.chromosome_bound_molecules_container import ChromosomeBoundMoleculeContainer

N_BASES = 1000
STRAND_MULTIPLICITY = 3

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		},
	'DNA polymerase':{
		},
	}

class Test_ChromosomeBoundMoleculeContainer(unittest.TestCase):
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
	def test_moleculesBound_all(self):
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
			self.container.moleculesBound()
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculesBound_byName(self):
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
			self.container.moleculesBound()
			)

		self.assertEqual(
			molecules,
			self.container.moleculesBound('DNA polymerase')
			)


	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_moleculesBound_position(self):
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
			self.container.moleculesBound(
				None,
				self.container.rootStrand(), 0, '+',
				200, 0
				)
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
			forks[0].attr('_sequencePosition'),
			startPosition
			)

		self.assertEqual(
			forks[1].attr('_sequencePosition'),
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
			forks[0].attr('_sequencePosition'),
			startPosition
			)

		self.assertEqual(
			forks[1].attr('_sequencePosition'),
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

		footprintParent = 4 # molecule origin is occupied by fork on parent
		regionParent = [111, 112, 113, 114]

		footprintChild = 2
		regionChild = [109, 110]

		self.container.moleculeLocationIsFork(mol, forkStop,
			forwardExtent, reverseExtent)

		chromosomeIndex = mol.attr('_globalIndex') + self.container._offset
		# TODO: make the above into a private method of the container class

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

		footprintParent = 4 # molecule origin is occupied by fork on parent
		regionParent = [96, 97, 98, 99]

		footprintChild = 2
		regionChild = [100, 101]

		self.container.moleculeLocationIsFork(mol, forkStart,
			forwardExtent, reverseExtent)

		chromosomeIndex = mol.attr('_globalIndex') + self.container._offset
		# TODO: make the above into a private method of the container class

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




def createContainer():
	container = ChromosomeBoundMoleculeContainer(N_BASES, STRAND_MULTIPLICITY,
		MOLECULE_ATTRIBUTES)

	return container
