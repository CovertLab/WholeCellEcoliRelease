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

import wholecell.utils.unique_objects_container

TEST_KB = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		}
	}


class Test_UniqueMoleculesContainer(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		self.container = wholecell.utils.unique_objects_container.UniqueObjectsContainer(
			TEST_KB)
		
		self.container.moleculesNew(
			'RNA polymerase',
			10,
			)

		self.container.moleculesNew(
			'RNA polymerase',
			5,
			boundToChromosome = True,
			chromosomeLocation = 0
			)

		self.container.moleculesNew(
			'RNA polymerase',
			5,
			boundToChromosome = True,
			chromosomeLocation = 50
			)


	def tearDown(self):
		pass

	@noseAttrib.attr('smalltest')
	def test_add_molecule(self):
		self.container.moleculeNew('RNA polymerase')

		self.assertEqual(len(self.container.molecules('RNA polymerase')), 21)


	@noseAttrib.attr('smalltest')
	def test_add_molecules(self):
		self.container.moleculesNew('RNA polymerase', 20)

		self.assertEqual(
			len(self.container.molecules('RNA polymerase')),
			40
			)


	@noseAttrib.attr('smalltest')
	def test_empty_query(self):
		molecules = self.container.evaluateQuery('RNA polymerase')

		self.assertEqual(len(molecules), 20)


	@noseAttrib.attr('smalltest')
	def test_bool_query(self):
		molecules = self.container.evaluateQuery(
			'RNA polymerase',
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
			'RNA polymerase',
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
			'RNA polymerase',
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
		for molecule in self.container.iterMolecules('RNA polymerase'):
			molecule.attrIs('boundToChromosome', True)
			molecule.attrIs('chromosomeLocation', 100)

		for molecule in self.container.iterMolecules('RNA polymerase'):
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
		query = self.container.queryNew('RNA polymerase', boundToChromosome = ('==', True))

		self.container.updateQueries()

		for molecule in query.iterMolecules():
			self.assertEqual(
				molecule.attr('boundToChromosome'),
				True
				)

		for molecule in query.iterMolecules():
			molecule.attrIs('boundToChromosome', False)

		self.container.updateQueries()

		self.assertEqual(query.molecules(), set())


	@noseAttrib.attr('smalltest')
	def test_delete_molecules(self):
		molecules = self.container.molecules('RNA polymerase')

		self.container.moleculesDel(molecules)

		self.assertEqual(
			self.container.molecules('RNA polymerase'),
			set()
			)

		# Make sure access to deleted molecules is blocked
		molecule = molecules.pop()

		with self.assertRaises(Exception) as context:
			molecule.attr('boundToChromosome')

		self.assertEqual(context.exception.message, 'Attempted to access an inactive molecule.')

		with self.assertRaises(Exception) as context:
			molecule.attrIs('boundToChromosome', False)

		self.assertEqual(context.exception.message, 'Attempted to access an inactive molecule.')


	@noseAttrib.attr('smalltest')
	def test_molecule_set_operations(self):
		allMolecules = self.container.molecules('RNA polymerase')

		chromosomeBound = self.container.evaluateQuery(
			'RNA polymerase',
			boundToChromosome = ('==', True)
			)

		self.assertEqual(allMolecules & allMolecules, allMolecules)

		self.assertEqual(allMolecules & chromosomeBound, chromosomeBound)

		self.assertEqual(allMolecules | chromosomeBound, allMolecules)

		self.assertEqual(len(allMolecules - chromosomeBound), 10)


	@noseAttrib.attr('smalltest')
	def test_time_setting(self):
		self.container._timeIs(50)

		allMolecules = self.container.molecules('RNA polymerase')
		newTime = self.container.evaluateQuery('RNA polymerase', _time = ('==', 50))

		self.assertEqual(allMolecules, newTime)


	@noseAttrib.attr('smalltest')
	def test_deleted_entry_flushing(self):
		# First, make sure that deleted entries are not overwritten
		molecules = self.container.molecules('RNA polymerase')
		indexes = {molecule._objectIndex for molecule in molecules}

		self.container.moleculesDel(molecules)

		newMolecule = self.container.moleculeNew(
			'RNA polymerase'
			)

		self.assertTrue(newMolecule._objectIndex not in indexes)

		# Next, flush the deleted entries and confirm that an old entry is overwritten
		self.container._flushDeleted()
		newMolecule = self.container.moleculeNew(
			'RNA polymerase'
			)

		self.assertTrue(newMolecule._objectIndex in indexes)


