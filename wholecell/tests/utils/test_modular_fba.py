#!/usr/bin/env python

"""
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/19/2013
"""

from __future__ import division

import unittest
import warnings

import numpy as np
import numpy.testing as npt
from wholecell.utils.modular_fba import FluxBalanceAnalysis

import nose.plugins.attrib as noseAttrib

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


	@noseAttrib.attr("smalltest", "fba")
	def test_instantiation(self):
		fba = FluxBalanceAnalysis(**_testStandard)


	@noseAttrib.attr("smalltest", "fba")
	def test_standard_IDs(self):
		fba = FluxBalanceAnalysis(**_testStandard)

		self.assertEqual(
			set(fba.externalMoleculeIDs()),
			{"A", "D"}
			)

		self.assertEqual(
			set(fba.outputMoleculeIDs()),
			{"B", "C", "E"}
			)


	@noseAttrib.attr("smalltest", "fba")
	def test_standard_noInput(self):
		fba = FluxBalanceAnalysis(**_testStandard)

		self.assertEqual(
			fba.biomassReactionFlux(),
			0
			)

		for moleculeID, change in zip(fba.outputMoleculeIDs(), fba.outputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(0, change)


	@noseAttrib.attr("smalltest", "fba")
	def test_standard(self):
		fba = FluxBalanceAnalysis(**_testStandard)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.externalMoleculeLevelsIs([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.externalMoleculeIDs()
			])

		self.assertEqual(
			fba.biomassReactionFlux(),
			1.0
			)

		for moleculeID, change in zip(fba.outputMoleculeIDs(), fba.outputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(10, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(10, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(20, change)


	@noseAttrib.attr("smalltest", "fba")
	def test_homeostatic_noInitial(self):
		fba = FluxBalanceAnalysis(**_testTargetMolecules)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.externalMoleculeLevelsIs([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.externalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":0,
			"C":0,
			"E":0,
			}

		fba.internalMoleculeLevelsIs([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.internalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.outputMoleculeIDs(), fba.outputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(10, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(10, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(20, change)


	@noseAttrib.attr("smalltest", "fba")
	def test_homeostatic_atObjective(self):
		fba = FluxBalanceAnalysis(**_testTargetMolecules)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.externalMoleculeLevelsIs([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.externalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":10,
			"C":10,
			"E":20,
			}

		fba.internalMoleculeLevelsIs([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.internalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.outputMoleculeIDs(), fba.outputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(0, change)


	@noseAttrib.attr("smalltest", "fba")
	def test_homeostatic_belowObjective(self):
		fba = FluxBalanceAnalysis(**_testTargetMolecules)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.externalMoleculeLevelsIs([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.externalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":5,
			"C":5,
			"E":10,
			}

		fba.internalMoleculeLevelsIs([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.internalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.outputMoleculeIDs(), fba.outputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(5, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(5, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(10, change)


	@noseAttrib.attr("smalltest", "fba")
	def test_homeostatic_singleBelowObjective(self):
		fba = FluxBalanceAnalysis(**_testTargetMolecules)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.externalMoleculeLevelsIs([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.externalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":5,
			"C":10,
			"E":20,
			}

		fba.internalMoleculeLevelsIs([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.internalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.outputMoleculeIDs(), fba.outputMoleculeLevelsChange()):
			if moleculeID == "B":
				self.assertAlmostEqual(5, change)

			elif moleculeID == "C":
				self.assertAlmostEqual(0, change)

			elif moleculeID == "E":
				self.assertAlmostEqual(0, change)


	@noseAttrib.attr("smalltest", "fba")
	def test_homeostatic_singleAboveObjective(self):
		fba = FluxBalanceAnalysis(
			**_testTargetMolecules
			)

		externalMoleculeLevels = {
			"A":50,
			"D":20
			}

		fba.externalMoleculeLevelsIs([
			externalMoleculeLevels[moleculeID]
			for moleculeID in fba.externalMoleculeIDs()
			])

		internalMoleculeLevels = {
			"B":15,
			"C":10,
			"E":20,
			}

		fba.internalMoleculeLevelsIs([
			internalMoleculeLevels[moleculeID]
			for moleculeID in fba.internalMoleculeIDs()
			])

		for moleculeID, change in zip(fba.outputMoleculeIDs(), fba.outputMoleculeLevelsChange()):
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
