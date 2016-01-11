#!/usr/bin/env python

"""
Test FlexTFbaModel.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/15/2013
"""

# TODO: DELETE THIS FILE
from __future__ import division

import unittest
import warnings

import numpy
# import wholecell.utils.flex_t_fba_model

import nose.plugins.attrib as noseAttrib


class Test_FlexTFbaModel(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		# Model comes from:
		# Covert MW, Schilling CH, Palsson BO. "Regulation of gene expression in flux balance models of metabolism" J Theor Biol. 2001. 213(1): 73-88.
		# [http://covertlab.stanford.edu/publicationpdfs/Covert2001b.pdf]
		# Modifications:
		# No external metabolites
		# R5a and R5b are identical (just different regulatory rules) --> merged to R5
		cls.metIds = [	"A[c]", "B[c]", "C[c]", "D[c]", "E[c]", "F[c]",
						"G[c]", "H[c]", "ATP[c]", "NADH[c]", "O2[c]"
					 ]
		cls.mediaEx = [
						{"rxnId": "Tc1", "met": "A[c]"},
						{"rxnId": "Tc2", "met": "A[c]"},
						{"rxnId": "Tf", "met": "F[c]"},
						{"rxnId": "Td", "met": "D[c]"},
						{"rxnId": "Te", "met": "E[c]"},
						{"rxnId": "Th", "met": "H[c]"},
						{"rxnId": "To2", "met": "O2[c]"},
					  ]
		cls.rxns = [
						{"id": "R1", "stoichiometry": [
							{"molecule": "A", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "B", "form": "m", "location": "c", "coeff": 1},
													  ]
						},
						{"id": "R2a", "stoichiometry": [
							{"molecule": "B", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "location": "c", "coeff": 2},
							{"molecule": "NADH", "form": "m", "location": "c", "coeff": 2},
							{"molecule": "C", "form": "m", "location": "c", "coeff": 1},
													  ]
						},
						{"id": "R2b", "stoichiometry": [
							{"molecule": "C", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "location": "c", "coeff": -2},
							{"molecule": "NADH", "form": "m", "location": "c", "coeff": -2},
							{"molecule": "B", "form": "m", "location": "c", "coeff": 1},
													  ]
						},
						{"id": "R3", "stoichiometry": [
							{"molecule": "B", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "F", "form": "m", "location": "c", "coeff": 1},
													  ]
						},
						{"id": "R4", "stoichiometry": [
							{"molecule": "C", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "G", "form": "m", "location": "c", "coeff": 1},
													  ]
						},
						{"id": "R5", "stoichiometry": [
							{"molecule": "G", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "C", "form": "m", "location": "c", "coeff": 0.8},
							{"molecule": "NADH", "form": "m", "location": "c", "coeff": 2},
													  ]
						},
						{"id": "R6", "stoichiometry": [
							{"molecule": "C", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "location": "c", "coeff": 2},
							{"molecule": "D", "form": "m", "location": "c", "coeff": 3},
													  ]
						},
						{"id": "R7", "stoichiometry": [
							{"molecule": "C", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "NADH", "form": "m", "location": "c", "coeff": -4},
							{"molecule": "E", "form": "m", "location": "c", "coeff": 3},
													  ]
						},
						{"id": "R8a", "stoichiometry": [
							{"molecule": "G", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "NADH", "form": "m", "location": "c", "coeff": -2},
							{"molecule": "H", "form": "m", "location": "c", "coeff": 1},
													  ]
						},
						{"id": "R8b", "stoichiometry": [
							{"molecule": "G", "form": "m", "location": "c", "coeff": 1},
							{"molecule": "ATP", "form": "m", "location": "c", "coeff": 1},
							{"molecule": "NADH", "form": "m", "location": "c", "coeff": 2},
							{"molecule": "H", "form": "m", "location": "c", "coeff": -1},
													  ]
						},
						{"id": "Rres", "stoichiometry": [
							{"molecule": "NADH", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "O2", "form": "m", "location": "c", "coeff": -1},
							{"molecule": "ATP", "form": "m", "location": "c", "coeff": 1},
													  ]
						}
				   ]
		cls.biomass = [
						{"id": "C[c]", "coeff": -1},
						{"id": "F[c]", "coeff": -1},
						{"id": "H[c]", "coeff": -1},
						{"id": "ATP[c]", "coeff": -10},
					  ]
		cls.atpId = "ATP[c]"
		cls.params = {"alpha": 10, "beta": 1000, "gamma": -1}
		

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		m = wholecell.utils.flex_t_fba_model.FlexTFbaModel(	metIds = self.metIds, rxns = self.rxns, mediaEx = self.mediaEx,
														biomass = self.biomass, atpId = self.atpId, params = self.params)
		
		self.lb = wholecell.utils.flex_t_fba_model.bounds(["thermodynamic", "exchange"], m.rxnIds(), False)
		self.ub = wholecell.utils.flex_t_fba_model.bounds(["thermodynamic", "exchange"], m.rxnIds(), True)
		
		# Exchange bounds given in paper
		self.lb.valuesIs(m.rxnIdxs(["mediaEx_Tc1"]), "exchange", -10.5)
		self.lb.valuesIs(m.rxnIdxs(["mediaEx_Tc2"]), "exchange", -10.5)
		self.lb.valuesIs(m.rxnIdxs(["mediaEx_Tf"]),  "exchange", -5.0)
		self.lb.valuesIs(m.rxnIdxs(["mediaEx_Th"]),  "exchange", -5.0)
		self.lb.valuesIs(m.rxnIdxs(["mediaEx_To2"]), "exchange", -15.0)

		self.ub.valuesIs(m.rxnIdxs(["mediaEx_Td"]), "exchange", 12.0)
		self.ub.valuesIs(m.rxnIdxs(["mediaEx_Te"]), "exchange", 12.0)

		# Reaction irreversibility
		self.lb.valuesIs(m.rxnGroup("real").idxs(), "thermodynamic", 0)

		# Biomass return is zero
		self.lb.valuesIs(m.rxnGroup("x").idxs(), "exchange", 0)

		m.v_lowerIs(idxs = m.rxnGroup("lowerMutable").idxs(), values = self.lb.mergedValues(m.rxnGroup("lowerMutable").idxs()))
		m.v_upperIs(idxs = m.rxnGroup("upperMutable").idxs(), values = self.ub.mergedValues(m.rxnGroup("upperMutable").idxs()))

		self.m = m

	def tearDown(self):
		pass

	# @noseAttrib.attr('smalltest')
	def test_wildType(self):
		m = self.m

		sol = m.solution()

		# Assert biomass matches expected (?) value (the value I saw)
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], 3.51818182, rtol = 0, atol = 1e-8))

		# Assert all f_i fluxes match g_bio (in the wild type case...)
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], sol[m.rxnGroup("f").idxs()], rtol = 0, atol = 1e-8))

		# Make sure the following has no runtime errors
		m.metaboliteProduction(m.metIdxs(m.metIds()), m.solution())

	# @noseAttrib.attr('smalltest')
	def test_KO_f(self):
		m = self.m 

		m.v_lowerIs(m.rxnIdxs(["rxn_R3"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R3"]), 0)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Tf"]), 0)
		m.v_upperIs(m.rxnIdxs(["mediaEx_Tf"]), 0)
		
		sol = m.solution()

		# F production should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_F[c]"])], 0, rtol = 0, atol = 1e-8))

		# Biomass should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], 0, rtol = 0, atol = 1e-8))

		# All non-F metabolites should be fine
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_C[c]", "f_H[c]", "f_ATP[c]"])], 3.51818182, rtol = 0, atol = 1e-8))

	# @noseAttrib.attr('smalltest')
	def test_KO_F(self):
		m = self.m 

		m.v_lowerIs(m.rxnIdxs(["rxn_R3"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R3"]), 0)
		m.v_lowerIs(m.rxnIdxs(["mediaEx_Tf"]), 0)
		m.v_upperIs(m.rxnIdxs(["mediaEx_Tf"]), 0)
		
		sol = m.solution()

		# F production should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_F[c]"])], 0, rtol = 0, atol = 1e-8))

		# Biomass should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], 0, rtol = 0, atol = 1e-8))

		# All non-F metabolites should be fine
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_C[c]", "f_H[c]", "f_ATP[c]"])], 3.51818182, rtol = 0, atol = 1e-8))

	# @noseAttrib.attr('smalltest')
	def test_KO_C(self):
		m = self.m 

		m.v_lowerIs(m.rxnIdxs(["rxn_R2a"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R2a"]), 0)
		m.v_lowerIs(m.rxnIdxs(["rxn_R5"]), 0)
		m.v_upperIs(m.rxnIdxs(["rxn_R5"]), 0)
		
		sol = m.solution()

		# Pretty lethal
		self.assertTrue(numpy.allclose(sol[m.rxnGroup("real").idxs()], 0, rtol = 0, atol = 1e-8))

	# @noseAttrib.attr('smalltest')
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
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_H[c]"])], 0, rtol = 0, atol = 1e-8))

		# Biomass should be zero
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["g_bio"])], 0, rtol = 0, atol = 1e-8))

		# All non-F metabolites should be fine
		self.assertTrue(numpy.allclose(sol[m.rxnIdxs(["f_C[c]", "f_F[c]", "f_ATP[c]"])], 3.51818182, rtol = 0, atol = 1e-8))