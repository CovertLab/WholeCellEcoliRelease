'''
test_unique_molecules_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 2/27/2014
'''

import unittest
import cPickle
import os

import numpy as np
import nose.plugins.attrib as noseAttrib

import wholecell.states.unique_molecules as wcUM

TEST_MOLECULE = 'RNA polymerase'

TEST_ATTRIBUTES = {
	'boundToChromosome':'bool',
	'chromosomeLocation':'uint32'
	}


class Test_UniqueMoleculesContainer(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		self.container = wcUM.UniqueMoleculesContainer(TEST_MOLECULE, TEST_ATTRIBUTES)
		
		self.container.moleculesNew(
			10,
			)

		self.container.moleculesNew(
			5,
			boundToChromosome = True,
			chromosomeLocation = 0
			)

		self.container.moleculesNew(
			5,
			boundToChromosome = True,
			chromosomeLocation = 50
			)


	def tearDown(self):
		pass


	@noseAttrib.attr('smalltest')
	def test_add_molecule(self):
		self.container.moleculeNew()

		self.assertEqual(len(self.container.molecules()), 21)


	@noseAttrib.attr('smalltest')
	def test_add_molecules(self):
		self.container.moleculesNew(20)

		self.assertEqual(
			len(self.container.molecules()),
			40
			)


	@noseAttrib.attr('smalltest')
	def test_empty_query(self):
		molecules = self.container.evaluateQuery()

		self.assertEqual(len(molecules), 20)


	@noseAttrib.attr('smalltest')
	def test_bool_query(self):
		molecules = self.container.evaluateQuery(
			boundToChromosome = ('==', False)
			)

		self.assertEqual(len(molecules), 10)

		for molecule in molecules:
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				False
				)


	@noseAttrib.attr('smalltest')
	def test_numeric_query(self):
		molecules = self.container.evaluateQuery(
			chromosomeLocation = ('>', 0)
			)

		self.assertEqual(len(molecules), 5)

		for molecule in molecules:
			self.assertGreater(
				molecule.attr('chromosomeLocation'),
				0
				)


	@noseAttrib.attr('smalltest')
	def test_compound_query(self):
		molecules = self.container.evaluateQuery(
			boundToChromosome = ('!=', False),
			chromosomeLocation = ('>', 0)
			)

		self.assertEqual(len(molecules), 5)

		for molecule in molecules:
			self.assertGreater(
				molecule.attr('chromosomeLocation'),
				0
				)


	@noseAttrib.attr('smalltest')
	def test_attribute_setting(self):
		for molecule in self.container.molecules():
			molecule.attrIs('boundToChromosome', True)
			molecule.attrIs('chromosomeLocation', 100)

		for molecule in self.container.molecules():
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				True
				)

			self.assertEqual(
				molecule.attr('chromosomeLocation'),
				100
				)


	@noseAttrib.attr('smalltest')
	def test_query_objects(self):
		query = self.container.queryNew(boundToChromosome = ('==', True))

		self.container.updateQueries()

		for molecule in query.molecules():
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				True
				)

		for molecule in query.molecules():
			molecule.attrIs('boundToChromosome', False)

		self.container.updateQueries()

		self.assertEqual(query.molecules(), set())


	@noseAttrib.attr('smalltest')
	def test_delete_molecules(self):
		molecules = self.container.molecules()

		self.container.moleculesDel(molecules)

		self.assertEqual(
			self.container.molecules(),
			set()
			)


