#!/usr/bin/env python

"""
Test flextFbaModel.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/15/2013
"""

import unittest
import warnings

import numpy
import wholecell.util.flextFbaModel

class Test_flextFbaModel(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		# Model comes from:
		# Covert MW, Schilling CH, Palsson BØ. "Regulation of gene expression in flux balance models of metabolism" J Theor Biol. 2001. 213(1): 73-88.
		# [http://covertlab.stanford.edu/publicationpdfs/Covert2001b.pdf]
		# Modifications:
		# No external metabolites
		# R5a and R5b are identical (just different regulatory rules) --> merged to R5
		cls.metIds = [	"A:m[c]", "B:m[c]", "C:m[c]", "D:m[c]", "E:m[c]", "F:m[c]",
						"G:m[c]", "H:m[c]", "ATP:m[c]", "NADH:m[c]", "O2:m[c]"
					 ]
		cls.mediaEx = [
						{"rxn": "Tc1", "met": "A:m[c]"},
						{"rxn": "Tc2", "met": "A:m[c]"},
						{"rxn": "Tf", "met": "F:m[c]"},
						{"rxn": "Td", "met": "D:m[c]"},
						{"rxn": "Te", "met": "E:m[c]"},
						{"rxn": "Th", "met": "H:m[c]"},
						{"rxn": "To2", "met": "O2:m[c]"},
					  ]
		cls.rxns = [
						{"id": "R1", "stoichiometry": [
							{"molecule": "A", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "B", "form": "m", "compartment": "c", "coeff": 1},
													  ]
						},
						{"id": "R2a", "stoichiometry": [
							{"molecule": "B", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "compartment": "c", "coeff": 2},
							{"molecule": "NADH", "form": "m", "compartment": "c", "coeff": 2},
							{"molecule": "C", "form": "m", "compartment": "c", "coeff": 1},
													  ]
						},
						{"id": "R2b", "stoichiometry": [
							{"molecule": "C", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "compartment": "c", "coeff": -2},
							{"molecule": "NADH", "form": "m", "compartment": "c", "coeff": -2},
							{"molecule": "B", "form": "m", "compartment": "c", "coeff": 1},
													  ]
						},
						{"id": "R3", "stoichiometry": [
							{"molecule": "B", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "F", "form": "m", "compartment": "c", "coeff": 1},
													  ]
						},
						{"id": "R4", "stoichiometry": [
							{"molecule": "C", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "G", "form": "m", "compartment": "c", "coeff": 1},
													  ]
						},
						{"id": "R5", "stoichiometry": [
							{"molecule": "G", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "C", "form": "m", "compartment": "c", "coeff": 0.8},
							{"molecule": "NADH", "form": "m", "compartment": "c", "coeff": 2},
													  ]
						},
						{"id": "R6", "stoichiometry": [
							{"molecule": "C", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "compartment": "c", "coeff": 2},
							{"molecule": "D", "form": "m", "compartment": "c", "coeff": 3},
													  ]
						},
						{"id": "R7", "stoichiometry": [
							{"molecule": "C", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "NADH", "form": "m", "compartment": "c", "coeff": -4},
							{"molecule": "E", "form": "m", "compartment": "c", "coeff": 3},
													  ]
						},
						{"id": "R8a", "stoichiometry": [
							{"molecule": "G", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "NADH", "form": "m", "compartment": "c", "coeff": -2},
							{"molecule": "H", "form": "m", "compartment": "c", "coeff": 1},
													  ]
						},
						{"id": "R8b", "stoichiometry": [
							{"molecule": "G", "form": "m", "compartment": "c", "coeff": 1},
							{"molecule": "ATP", "form": "m", "compartment": "c", "coeff": 1},
							{"molecule": "NADH", "form": "m", "compartment": "c", "coeff": 2},
							{"molecule": "H", "form": "m", "compartment": "c", "coeff": -1},
													  ]
						},
						{"id": "Rres", "stoichiometry": [
							{"molecule": "NADH", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "O2", "form": "m", "compartment": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "compartment": "c", "coeff": 1},
													  ]
						}
				   ]
		cls.biomass = [
						{"id": "C:m[c]", "r": -1},
						{"id": "F:m[c]", "r": -1},
						{"id": "H:m[c]", "r": -1},
						{"id": "ATP:m[c]", "r": -10},
					  ]
		cls.atpId = "ATP:m[c]"
		cls.params = {"alpha": 10, "beta": 1000, "gamma": -1}

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		m = wholecell.util.flextFbaModel.flextFbaModel(	metIds = self.metIds, rxns = self.rxns, mediaEx = self.mediaEx,
														biomass = self.biomass, atpId = self.atpId, params = self.params)
		
		# Given in paper
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Tc1"]), -10.5)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Tc2"]), -10.5)
		m.v_upperIs(m.rxnIdxs(["mediaEx_Td"]),   12.0)
		m.v_upperIs(m.rxnIdxs(["mediaEx_Te"]),   12.0)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Tf"]),  -5.00)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Th"]),  -5.00)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_To2"]), -15.0)

		# Irreversibility
		m.v_lowerIs(m.rxnGroup("real").idxs(), 0)			# Irreversibility

		# Biomass return is zero
		m.v_lowerIs(m.rxnGroup("x").idxs(), 0)				# Biomass doesn't return anything

		self.m = m

	def tearDown(self):
		pass

	def test_wildType(self):
		m = self.m

		sol = m.solution()

		# Assert biomass matches expected (?) value (the value I saw)
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], 3.51818182, rtol = 0, atol = 1e-8))

		# Assert all f_i fluxes match g_bio (in the wild type case...)
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], sol[m.rxnGroup("f").idxs()], rtol = 0, atol = 1e-8))

	def test_KO_f(self):
		m = self.m 

		m.v_lowerIs(m.rxnIdxs(["rxn_R3"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R3"]), 0)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Tf"]), 0)
		m.v_upperIs(m.rxnIdxs(["mediaEx_Tf"]), 0)
		
		sol = m.solution()

		# F production should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_F:m[c]"])], 0, rtol = 0, atol = 1e-8))

		# Biomass should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], 0, rtol = 0, atol = 1e-8))

		# All non-F metabolites should be fine
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_C:m[c]", "f_H:m[c]", "f_ATP:m[c]"])], 3.51818182, rtol = 0, atol = 1e-8))

	def test_KO_F(self):
		m = self.m 

		m.v_lowerIs(m.rxnIdxs(["rxn_R3"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R3"]), 0)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Tf"]), 0)
		m.v_upperIs(m.rxnIdxs(["mediaEx_Tf"]), 0)
		
		sol = m.solution()

		# F production should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_F:m[c]"])], 0, rtol = 0, atol = 1e-8))

		# Biomass should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], 0, rtol = 0, atol = 1e-8))

		# All non-F metabolites should be fine
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_C:m[c]", "f_H:m[c]", "f_ATP:m[c]"])], 3.51818182, rtol = 0, atol = 1e-8))

	def test_KO_C(self):
		m = self.m 

		m.v_lowerIs(m.rxnIdxs(["rxn_R2a"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R2a"]), 0)
		m.v_lowerIs(m.rxnIdxs(["rxn_R5"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R5"]), 0)
		
		sol = m.solution()

		# Pretty lethal
		self.assertTrue(numpy.allclose(sol[m.rxnGroup("real").idxs()], 0, rtol = 0, atol = 1e-8))

	def test_KO_H(self):
		m = self.m 

		m.v_lowerIs(m.rxnIdxs(["rxn_R8a"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R8a"]), 0)
		m.v_lowerIs(m.rxnIdxs(["rxn_R8b"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R8b"]), 0)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Th"]), 0)
		m.v_upperIs(m.rxnIdxs(["mediaEx_Th"]), 0)
		
		sol = m.solution()

		# H production should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_H:m[c]"])], 0, rtol = 0, atol = 1e-8))

		# Biomass should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], 0, rtol = 0, atol = 1e-8))

		# All non-F metabolites should be fine
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_C:m[c]", "f_F:m[c]", "f_ATP:m[c]"])], 3.51818182, rtol = 0, atol = 1e-8))