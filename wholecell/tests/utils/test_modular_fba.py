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

		fba.run()

		self.assertEqual(
			fba.objectiveReactionFlux(),
			0
			)

		self.assertEqual(
			fba.outputMoleculeLevelsChange().tolist(),
			[0, 0, 0]
			)


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

		fba.run()

		self.assertEqual(
			fba.objectiveReactionFlux(),
			1.0
			)

		npt.assert_allclose(
			fba.outputMoleculeLevelsChange(),
			[10, 10, 20]
			)

# TODO: tests for enzymes
# TODO: tests for mass accumulation
# TODO: tests for flexible FBA
# TODO: tests for exceptions
