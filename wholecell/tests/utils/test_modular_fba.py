#!/usr/bin/env python

"""
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/19/2013
"""

from __future__ import absolute_import, division, print_function

import unittest

from wholecell.utils.modular_fba import FluxBalanceAnalysis
from six.moves import zip


# TODO: test all solvers

_testStandard = dict(
	reactionStoich = {
		"A to B":{"A":-1, "B":+1},
		"2B to C":{"B":-2, "C":+1},
		"A + D to E":{"A":-1, "D":-1, "E":+1},
		},
	externalExchangedMolecules = ["A", "D"],
	objective = {
		"B":10,
		"C":10,
		"E":20,
		}
	)

_testTargetMolecules = dict(
	reactionStoich = {
		"A to B":{"A":-1, "B":+1},
		"B to A":{"A":+1, "B":-1},
		"2B to C":{"B":-2, "C":+1},
		"C to 2B":{"B":+2, "C":-1},
		"A + D to E":{"A":-1, "D":-1, "E":+1},
		"E to A + D":{"A":+1, "D":+1, "E":-1},
		},
	externalExchangedMolecules = ["A", "D"],
	objective = {
		"B":10,
		"C":10,
		"E":20,
		},
	objectiveType = "homeostatic"
	)

class Test_FluxBalanceAnalysis(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		pass


	def tearDown(self):
		pass


	def test_instantiation(self):
		fba = FluxBalanceAnalysis(**_testStandard)


	def test_standard_IDs(self):
		fba = FluxBalanceAnalysis(**_testStandard)

		self.assertEqual(
			set(fba.getExternalMoleculeIDs()),
			{"A", "D"}
			)

		self.assertEqual(
			set(fba.getOutputMoleculeIDs()),
			{"B", "C", "E"}
			)


	def test_standard_noInput(self):
		fba = FluxBalanceAnalysis(**_testStandard)

		self.assertEqual(
			fba.getBiomassReactionFlux(),
			0
			)

		for moleculeID, change in zip(fba.getOutputMoleculeIDs(), fba.getOutputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(0, change)


	def test_standard(self):
		fba = FluxBalanceAnalysis(**_testStandard)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.setExternalMoleculeLevels([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.getExternalMoleculeIDs()
			])

		self.assertEqual(
			fba.getBiomassReactionFlux(),
			1.0
			)

		for moleculeID, change in zip(fba.getOutputMoleculeIDs(), fba.getOutputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(10, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(10, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(20, change)


	def test_homeostatic_noInitial(self):
		fba = FluxBalanceAnalysis(**_testTargetMolecules)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.setExternalMoleculeLevels([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.getExternalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":0,
			"C":0,
			"E":0,
			}

		fba.setInternalMoleculeLevels([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.getInternalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.getOutputMoleculeIDs(), fba.getOutputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(10, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(10, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(20, change)


	def test_homeostatic_atObjective(self):
		fba = FluxBalanceAnalysis(**_testTargetMolecules)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.setExternalMoleculeLevels([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.getExternalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":10,
			"C":10,
			"E":20,
			}

		fba.setInternalMoleculeLevels([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.getInternalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.getOutputMoleculeIDs(), fba.getOutputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(0, change)


	def test_homeostatic_belowObjective(self):
		fba = FluxBalanceAnalysis(**_testTargetMolecules)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.setExternalMoleculeLevels([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.getExternalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":5,
			"C":5,
			"E":10,
			}

		fba.setInternalMoleculeLevels([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.getInternalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.getOutputMoleculeIDs(), fba.getOutputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(5, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(5, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(10, change)


	def test_homeostatic_singleBelowObjective(self):
		fba = FluxBalanceAnalysis(**_testTargetMolecules)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.setExternalMoleculeLevels([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.getExternalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":5,
			"C":10,
			"E":20,
			}

		fba.setInternalMoleculeLevels([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.getInternalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.getOutputMoleculeIDs(), fba.getOutputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(5, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(0, change)


	def test_homeostatic_singleAboveObjective(self):
		fba = FluxBalanceAnalysis(
			**_testTargetMolecules
			)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.setExternalMoleculeLevels([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.getExternalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":15,
			"C":10,
			"E":20,
			}

		fba.setInternalMoleculeLevels([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.getInternalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.getOutputMoleculeIDs(), fba.getOutputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(-5, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(0, change)

# TODO: tests for enzymes
# TODO: tests for mass accumulation
# TODO: tests for flexible FBA
# TODO: tests for exceptions
