"""
Test KnowledgeBase.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

import unittest
import warnings

import numpy
import wholecell.kb.KnowledgeBase

class Test_randStream(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		# TODO: Load fixture instead
		self.kb = wholecell.kb.KnowledgeBase.KnowledgeBase(dataFileName = "data/KnowledgeBase.xlsx",
															 seqFileName = "data/KnowledgeBase.fna")

	def tearDown(self):
		pass

	def test_construction(self):
		wholecell.kb.KnowledgeBase.KnowledgeBase(dataFileName = "data/KnowledgeBase.xlsx",
													seqFileName = "data/KnowledgeBase.fna")

	def test_metabolites(self):
		kb = self.kb

		self.assertEqual(722, len(kb.metabolites))

		self.assertEqual(82, len([x for x in kb.metabolites if x["hydrophobic"] == True]))
		self.assertEqual(640, len([x for x in kb.metabolites if x["hydrophobic"] == False]))

		met = next((x for x in kb.metabolites if x["id"] == "ACCOA"), None)
		self.assertNotEqual(met, None)
		self.assertEqual(dict, type(met))
		self.assertEqual("ACCOA", met["id"])
		self.assertEqual("Acetyl-CoA", met["name"])
		self.assertEqual("C23H34N7O17P3S1", met["formula"])
		self.assertEqual("CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N", met["smiles"])
		self.assertEqual(-4, met["charge"])
		self.assertEqual(805.538, met["mw"])
		self.assertEqual(False, met["hydrophobic"])
		self.assertEqual(0, met["mediaConc"])
		self.assertAlmostEqual(3.5812e+03, met["biomassConc"], places = 1)
		self.assertAlmostEqual(3.5812e+03, met["metabolismNewFlux"], places = 1)
		self.assertEqual(0, met["metabolismRecyclingFlux"])

		met = next((x for x in kb.metabolites if x["id"] == "AC"), None)
		self.assertAlmostEqual(0.304758, met["mediaConc"], places = 6)
		
	def test_geneticCode(self):
		kb = self.kb

		self.assertEqual(4, kb.translationTable)

	def test_genome(self):
		kb = self.kb

		self.assertEqual(580076, len(kb.genomeSeq))
		self.assertEqual("ACGT", "".join(sorted(list(set(kb.genomeSeq)))))