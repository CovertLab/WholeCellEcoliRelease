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

	#@noseAttrib.attr('focusTest')
	def test_parseReactionScript(self):
		rp = r.reactionParser()



		# Simple OR
		line = '( b0241  or  b0929  or  b1377  or  b2215 )'
		self.assertEqual(['MONOMER0-282', 'EG10671-MONOMER', 'G6700-MONOMER', 'EG10670-MONOMER'], rp.findEnzyme(line, 'hi')['enzymes'])

		# Simple and/OR
		line = '( ( b3670  and  b3671 )  or  ( b0077  and  b0078 ) )'
		self.assertEqual(['ACETOLACTSYNI-CPLX', 'ACETOLACTSYNIII-CPLX'], rp.findEnzyme(line, 'hi')['enzymes'])

		# Complicated complex with OR
		line = '( ( ( b3736  and  b3737  and  b3738 )  and  ( b3731  and  b3732  and  b3733  and  b3734  and  b3735 ) )  or  ( ( b3736  and  b3737  and  b3738 )  and  ( b3731  and  b3732  and  b3733  and  b3734  and  b3735 )  and  b3739 ) )'
		enzymes = rp.findEnzyme(line, 'hi')['enzymes']
		self.assertEqual(['ATPSYN-CPLX','UNKNOWN'], enzymes)

	def test_findEnzyme(self):
		rp = r.reactionParser()

		# Basic reaction with 'or'
		line = '( b0241  or  b0929  or  b1377  or  b2215 )'
		self.assertEqual([['MONOMER0-282'], 'or', ['EG10671-MONOMER'], 'or', ['G6700-MONOMER'], 'or', ['EG10670-MONOMER']], rp.findEnzyme(line)['enzymes'])

		# Basic complex with and
		line = '( b1252  and  b3005  and  b3006 )'
		self.assertEqual([['CPLX0-1923']], rp.findEnzyme(line)['enzymes'])
		line = '( b2677  and  b2678  and  b2679 )'
		self.assertEqual([['ABC-26-CPLX']], rp.findEnzyme(line)['enzymes'])

	@noseAttrib.attr('focusTest')
	def test_findEnzymeManual(self):
		rp = r.reactionParser()

		line = 'ATPASE-1-CPLX'

		line = '((PYRUVFORMLY-CPLX and EG11784-MONOMER) or PYRUVFORMLY-CPLX or KETOBUTFORMLY-MONOMER)'

	def test_iterateTree(self):
		rp = r.reactionParser()

		cmplxFrameId = 'C1'
		monomers = []
		monomerOrComplexToComplex = {'C1' : ['M1', 'M2']}
		rp.iterateTree(cmplxFrameId, monomers, monomerOrComplexToComplex)
		self.assertEqual(['M1', 'M2'], monomers)

		cmplxFrameId = 'C1'
		monomers = []
		monomerOrComplexToComplex = {'C1' : ['M1', 'C2'], 'C2' : ['M2', 'M3']}
		rp.iterateTree(cmplxFrameId, monomers, monomerOrComplexToComplex)
		self.assertEqual(['M1','M2','M3'].sort(), monomers.sort())


		cmplxFrameId = 'C1'
		monomers = []
		monomerOrComplexToComplex = {'C1' : ['M1', 'C2'], 'C2' : ['C3', 'C4'], 'C3' : ['M2','M3'], 'C4' : ['M4','M5']}
		rp.iterateTree(cmplxFrameId, monomers, monomerOrComplexToComplex)
		self.assertEqual(['M1','M2','M3','M4','M5'].sort(), monomers.sort())
