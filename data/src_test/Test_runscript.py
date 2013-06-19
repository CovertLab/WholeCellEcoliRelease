#!/usr/bin/env python
import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest

import data.src.runscript as r

import ipdb

class Test_Simulation(unittest.TestCase):

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

	def test_runscript(self):
		r.main()
	
	def test_getminCoord(self):
		gL = []
		for i in range(3):
			gL.append(r.gene())
			gL[i].left = i*10 + 10
		minCoord = r.getMinCoord(geneList = gL)
		self.assertEqual(minCoord, 10)

	def tet_getmaxCoord(self):
		gL = []
		for i in range(3):
			gL.append(r.gene())
			gL[i].left = i*10 + 10
		maxCoord = r.getMinCoord(geneList = gL)
		self.assertEqual(maxCoord, 30)

	def test_parseSigmaFactors(self):
		self.assertEqual(['S','D'], r.parseSigmaFactors('(RNA polymerase, sigma S (sigma 38) factor RNA polymerase, sigma 70 (sigma D) factor)'))

	def test_calculateWeight(self):
		met = r.metabolite()
		self.assertLess(abs(180.16 - met.calculateWeight('C6H12O6')), 0.006)
		self.assertLess(abs(18.01528 - met.calculateWeight('H2O')), 0.006)

	@noseAttrib.attr('focusTest')
	def test_parseReactionScript(self):
		rp = r.reactionParser()

		# Basic complex with and
		line = '( b1252  and  b3005  and  b3006 )'
		self.assertEqual(['CPLX0-1923'], rp.parseBracket(line))
		line = '( b2677  and  b2678  and  b2679 )'
		self.assertEqual(['ABC-26-CPLX'], rp.parseBracket(line))

		# Basic or
		line = '( b0241  or  b0929  or  b1377  or  b2215 )'
		self.assertEqual(['EG10670-MONOMER', 'EG10671-MONOMER', 'G6700-MONOMER', 'MONOMER0-282'], rp.parseBracket(line))







		line = '( ( b3670  and  b3671 )  or  ( b0077  and  b0078 ) )'

		nick = rp.parseRecursiveBracket(line[1:-1], [])
		ipdb.set_trace()


