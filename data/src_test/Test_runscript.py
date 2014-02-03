#!/usr/bin/env python
import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest
import os
import csv

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

	@noseAttrib.attr('parseTest')
	def test_runscript(self):
		r.main()

	@noseAttrib.attr('parseTest')
	def test_getminCoord(self):
		gL = []
		for i in range(3):
			gL.append(r.gene())
			gL[i].left = i*10 + 10
		minCoord = r.getMinCoord(geneList = gL)
		self.assertEqual(minCoord, 10)

	@noseAttrib.attr('parseTest')
	def tet_getmaxCoord(self):
		gL = []
		for i in range(3):
			gL.append(r.gene())
			gL[i].left = i*10 + 10
		maxCoord = r.getMinCoord(geneList = gL)
		self.assertEqual(maxCoord, 30)

	@noseAttrib.attr('parseTest')
	def test_parseSigmaFactors(self):
		self.assertEqual(['S','D'], r.parseSigmaFactors('(RNA polymerase, sigma S (sigma 38) factor RNA polymerase, sigma 70 (sigma D) factor)'))

	@noseAttrib.attr('parseTest')
	def test_calculateWeight(self):
		met = r.metabolite()
		self.assertLess(abs(180.16 - met.calculateWeight('C6H12O6')), 0.006)
		self.assertLess(abs(18.01528 - met.calculateWeight('H2O')), 0.006)

	@noseAttrib.attr('parseTest')
	def test_allComplexesCreated(self):
		Ecocyc_complexFrameIds = []
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_complexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			for row in csvreader:
				Ecocyc_complexFrameIds.append(row[0])
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_rna_protein_complexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			for row in csvreader:
				Ecocyc_complexFrameIds.append(row[0])
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw', 'Ecocyc_protein_small_molecule_complexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			for row in csvreader:
				Ecocyc_complexFrameIds.append(row[0])
		Ecocyc_complexFrameIds = set(Ecocyc_complexFrameIds)

		Parsed_complexFrameIds = []
		with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'),'rb') as csvfile:
			csvreader = csv.reader(csvfile, delimiter='\t', quotechar='"')
			csvreader.next()
			for row in csvreader:
				Parsed_complexFrameIds.append(row[0])
		Parsed_complexFrameIds = set(Parsed_complexFrameIds)
		self.assertEqual(len(Parsed_complexFrameIds.difference(Ecocyc_complexFrameIds)), 0)

	@noseAttrib.attr('parseTest')
	def test_findEnzyme(self):
		rp = r.reactionParser()

		# Basic reaction with 'or'
		line = '( b0241  or  b0929  or  b1377  or  b2215 )'
		self.assertEqual(['MONOMER0-282','EG10671-MONOMER','G6700-MONOMER','EG10670-MONOMER'], rp.findEnzyme(line)['enzymes'])

		# Basic complex with and
		line = '( b1252  and  b3005  and  b3006 )'
		self.assertEqual([['CPLX0-1923']], rp.findEnzyme(line)['enzymes'])
		line = '( b2677  and  b2678  and  b2679 )'
		self.assertEqual([['ABC-26-CPLX']], rp.findEnzyme(line)['enzymes'])

	@noseAttrib.attr('parseTest')
	def test_findEnzymeManual(self):
		rp = r.reactionParser()

		line = 'ATPASE-1-CPLX'
		self.assertEqual(['ATPASE-1-CPLX'], rp.findEnzymeManualCuration(line))

		line = '((FHLMULTI-CPLX and G7307-MONOMER) or (CPLX0-250 and G7307-MONOMER))'
		self.assertEqual([['FHLMULTI-CPLX', 'G7307-MONOMER'], ['CPLX0-250', 'G7307-MONOMER']],rp.findEnzymeManualCuration(line))

		line = '((PYRUVFORMLY-CPLX and EG11784-MONOMER) or PYRUVFORMLY-CPLX or KETOBUTFORMLY-MONOMER)'
		self.assertEqual([['PYRUVFORMLY-CPLX', 'EG11784-MONOMER'], ['PYRUVFORMLY-CPLX'], ['KETOBUTFORMLY-MONOMER']], rp.findEnzymeManualCuration(line))

		line = '((SAPD-MONOMER and TRKA-MONOMER and TRKG-MONOMER) or (SAPD-MONOMER and TRKA-MONOMER and TRKH-MONOMER))'
		self.assertEqual([['SAPD-MONOMER', 'TRKA-MONOMER', 'TRKG-MONOMER'], ['SAPD-MONOMER', 'TRKA-MONOMER', 'TRKH-MONOMER']], rp.findEnzymeManualCuration(line))

	@noseAttrib.attr('parseTest')
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

	@noseAttrib.attr('parseTest')
	def test_getEcocycComplexComponents(self):
		complexName = 'CPLX0-221'

		info = r.getEcocycComplexComponents(complexName)
		self.assertEqual(len(info), 2)
		self.assertEqual(info[0], ('APORNAP-CPLX', '1'))
		self.assertEqual(info[1], ('PD00440', '1'))

	@noseAttrib.attr('parseTest')
	def test_buildReactionInstanceFromClassList(self):
		# Reaction looks like:
		# tRNAval + L-valine + ATP + H+ -> L-valyl-tRNAval + AMP + diphosphate

		reaction = r.getEcocycReactionStoich('VALINE--TRNA-LIGASE-RXN')
		unmodified_form = 'RNA0-300'
		modified_form = 'RNA0-311'
		components_children = r.buildReactionInstanceFromClassList(reaction, modified_form, unmodified_form)
		components_children_test = [[{'classid' : 'VAL-tRNAs', 'instanceid' : 'RNA0-300'}],
									[{'classid' : 'Charged-VAL-tRNAs', 'instanceid' : 'RNA0-311'}]]

		components_children.sort()
		components_children[0].sort()
		components_children[1].sort()
		components_children_test.sort()
		components_children_test[0].sort()
		components_children_test[1].sort()


		self.assertEqual(components_children, components_children_test)

	@noseAttrib.attr('parseTest')
	def test_buildInstanceReaction(self):
		cart_product = ({'classid': 'VAL-tRNAs', 'instanceid': 'RNA0-300'}, {'classid': 'Charged-VAL-tRNAs', 'instanceid': 'RNA0-311'})
		reaction = r.getEcocycReactionStoich('VALINE--TRNA-LIGASE-RXN')

		new_rxn = r.buildInstanceReaction(cart_product, reaction)

		new_rxn_test = {'direction': 'LEFT-TO-RIGHT',
						'enzyme': ['VALS-MONOMER'],
						'components': [{'coeff': '-1', 'id': 'PROTON', 'isclass': False},
										{'coeff': '-1', 'id': 'RNA0-300', 'isclass': False},
										{'coeff': '-1', 'id': 'VAL', 'isclass': False},
										{'coeff': '-1', 'id': 'ATP', 'isclass': False},
										{'coeff': '1', 'id': 'RNA0-311', 'isclass': False},
										{'coeff': '1', 'id': 'PPI', 'isclass': False},
										{'coeff': '1', 'id': 'AMP', 'isclass': False}],
						'id': 'VALINE--TRNA-LIGASE-RXN'}

		self.assertEqual(new_rxn['direction'], new_rxn_test['direction'])
		self.assertEqual(new_rxn['enzyme'], new_rxn_test['enzyme'])
		self.assertEqual(new_rxn['id'], new_rxn_test['id'])
		self.assertEqual(new_rxn['components'].sort(), new_rxn_test['components'].sort())

	@noseAttrib.attr('parseTest')
	def test_getFormationReactions(self):
		unmodified_form = 'RNA0-300'
		frameId = 'RNA0-311'
		formation_reactions = r.getFormationReactions(frameId, unmodified_form)
		
		self.assertEqual(len(formation_reactions), 1)
		new_rxn_test = {'direction': 'LEFT-TO-RIGHT',
						'enzyme': ['VALS-MONOMER'],
						'components': [{'coeff': '-1', 'id': 'PROTON', 'isclass': False},
										{'coeff': '-1', 'id': 'RNA0-300', 'isclass': False},
										{'coeff': '-1', 'id': 'VAL', 'isclass': False},
										{'coeff': '-1', 'id': 'ATP', 'isclass': False},
										{'coeff': '1', 'id': 'RNA0-311', 'isclass': False},
										{'coeff': '1', 'id': 'PPI', 'isclass': False},
										{'coeff': '1', 'id': 'AMP', 'isclass': False}],
						'id': 'VALINE--TRNA-LIGASE-RXN'}
		self.assertEqual(formation_reactions[0]['direction'], new_rxn_test['direction'])
		self.assertEqual(formation_reactions[0]['enzyme'], new_rxn_test['enzyme'])
		self.assertEqual(formation_reactions[0]['id'], new_rxn_test['id'])
		self.assertEqual(formation_reactions[0]['components'].sort(), new_rxn_test['components'].sort())

