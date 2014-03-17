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

from wholecell.utils.chromosome_bound_molecules_container import ChromosomeBoundMoleculesContainer

N_BASES = 1000
STRAND_MULTIPLICITY = 3

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		},
	'DNA polymerase':{
		},
	}


class Test_UniqueObjectsContainer(unittest.TestCase):
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

	# Adding/removing objects
	@noseAttrib.attr('smalltest', 'chromosome', 'containerObject')
	def test_bind_molecule(self):
		mol = self.container.moleculeNew('DNA polymerase')

		self.container.moleculeLocation(mol, self.container.rootStrand(), 2000,
			'+', 50, 10)


def createContainer():
	container = ChromosomeBoundMoleculesContainer(N_BASES, STRAND_MULTIPLICITY,
		MOLECULE_ATTRIBUTES)

	return container
