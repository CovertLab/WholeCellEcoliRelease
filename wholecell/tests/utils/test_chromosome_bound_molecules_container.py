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


def createContainer():
	container = ChromosomeBoundMoleculeContainer(N_BASES, STRAND_MULTIPLICITY,
		MOLECULE_ATTRIBUTES)

	return container
