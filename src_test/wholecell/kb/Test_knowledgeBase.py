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

	def test_genes(self):
		kb = self.kb

		self.assertEqual(525, len(kb.genes))
		self.assertEqual(482, len([x for x in kb.genes if x["type"] == "mRNA"]))
		self.assertEqual(3, len([x for x in kb.genes if x["type"] == "rRNA"]))
		self.assertEqual(4, len([x for x in kb.genes if x["type"] == "sRNA"]))
		self.assertEqual(36, len([x for x in kb.genes if x["type"] == "tRNA"]))

		gene = next((x for x in kb.genes if x["id"] == "MG_001"), None)
		self.assertNotEqual(gene, None)
		self.assertEqual(dict, type(gene))
		self.assertEqual("MG_001", gene["id"])
		self.assertEqual("DNA polymerase III, beta subunit", gene["name"])
		self.assertEqual("dnaN", gene["symbol"])
		self.assertEqual("mRNA", gene["type"])
		self.assertEqual(686, gene["start"])
		self.assertEqual(1143, gene["len"])
		self.assertTrue(gene["dir"])
		self.assertEqual(
                "ATGAAAATATTAATTAATAAAAGTGAATTGAATAAAATTTTGAAAAAAAT" +
                "GAATAACGTTATTATTTCCAATAACAAAATAAAACCACATCATTCATATT" +
                "TTTTAATAGAGGCAAAAGAAAAAGAAATAAACTTTTATGCTAACAATGAA" +
                "TACTTTTCTGTCAAATGTAATTTAAATAAAAATATTGATATTCTTGAACA" +
                "AGGCTCCTTAATTGTTAAAGGAAAAATTTTTAACGATCTTATTAATGGCA" +
                "TAAAAGAAGAGATTATTACTATTCAAGAAAAAGATCAAACACTTTTGGTT" +
                "AAAACAAAAAAAACAAGTATTAATTTAAACACAATTAATGTGAATGAATT" +
                "TCCAAGAATAAGGTTTAATGAAAAAAACGATTTAAGTGAATTTAATCAAT" +
                "TCAAAATAAATTATTCACTTTTAGTAAAAGGCATTAAAAAAATTTTTCAC" +
                "TCAGTTTCAAATAATCGTGAAATATCTTCTAAATTTAATGGAGTAAATTT" +
                "CAATGGATCCAATGGAAAAGAAATATTTTTAGAAGCTTCTGACACTTATA" +
                "AACTATCTGTTTTTGAGATAAAGCAAGAAACAGAACCATTTGATTTCATT" +
                "TTGGAGAGTAATTTACTTAGTTTCATTAATTCTTTTAATCCTGAAGAAGA" +
                "TAAATCTATTGTTTTTTATTACAGAAAAGATAATAAAGATAGCTTTAGTA" +
                "CAGAAATGTTGATTTCAATGGATAACTTTATGATTAGTTACACATCGGTT" +
                "AATGAAAAATTTCCAGAGGTAAACTACTTTTTTGAATTTGAACCTGAAAC" +
                "TAAAATAGTTGTTCAAAAAAATGAATTAAAAGATGCACTTCAAAGAATTC" +
                "AAACTTTGGCTCAAAATGAAAGAACTTTTTTATGCGATATGCAAATTAAC" +
                "AGTTCTGAATTAAAAATAAGAGCTATTGTTAATAATATCGGAAATTCTCT" +
                "TGAGGAAATTTCTTGTCTTAAATTTGAAGGTTATAAACTTAATATTTCTT" +
                "TTAACCCAAGTTCTCTATTAGATCACATAGAGTCTTTTGAATCAAATGAA" +
                "ATAAATTTTGATTTCCAAGGAAATAGTAAGTATTTTTTGATAACCTCTAA" +
                "AAGTGAACCTGAACTTAAGCAAATATTGGTTCCTTCAAGATAA",
                gene["seq"]
			)
		self.assertEqual("MG_001", gene["rnaId"])

	def test_rnas(self):
		kb = self.kb

		self.assertEqual(525, len(kb.rnas))

		rna = next((x for x in kb.rnas if x["id"] == "MG_001"), None)
		self.assertNotEqual(rna, None)
		self.assertEqual(dict, type(rna))
		self.assertEqual("MG_001", rna["id"])
		self.assertEqual("DNA polymerase III, beta subunit", rna["name"])
		self.assertEqual("mRNA", rna["type"])
		self.assertAlmostEqual(8.8983e-5, rna["exp"], places = 9)
		self.assertAlmostEqual(146.9388, rna["halfLife"], places = 4)
		self.assertEqual(
			    "UACUUUUAUAAUUAAUUAUUUUCACUUAACUUAUUUUAAAACUUUUUUUA" +
                "CUUAUUGCAAUAAUAAAGGUUAUUGUUUUAUUUUGGUGUAGUAAGUAUAA" +
                "AAAAUUAUCUCCGUUUUCUUUUUCUUUAUUUGAAAAUACGAUUGUUACUU" +
                "AUGAAAAGACAGUUUACAUUAAAUUUAUUUUUAUAACUAUAAGAACUUGU" +
                "UCCGAGGAAUUAACAAUUUCCUUUUUAAAAAUUGCUAGAAUAAUUACCGU" +
                "AUUUUCUUCUCUAAUAAUGAUAAGUUCUUUUUCUAGUUUGUGAAAACCAA" +
                "UUUUGUUUUUUUUGUUCAUAAUUAAAUUUGUGUUAAUUACACUUACUUAA" +
                "AGGUUCUUAUUCCAAAUUACUUUUUUUGCUAAAUUCACUUAAAUUAGUUA" +
                "AGUUUUAUUUAAUAAGUGAAAAUCAUUUUCCGUAAUUUUUUUAAAAAGUG" +
                "AGUCAAAGUUUAUUAGCACUUUAUAGAAGAUUUAAAUUACCUCAUUUAAA" +
                "GUUACCUAGGUUACCUUUUCUUUAUAAAAAUCUUCGAAGACUGUGAAUAU" +
                "UUGAUAGACAAAAACUCUAUUUCGUUCUUUGUCUUGGUAAACUAAAGUAA" +
                "AACCUCUCAUUAAAUGAAUCAAAGUAAUUAAGAAAAUUAGGACUUCUUCU" +
                "AUUUAGAUAACAAAAAAUAAUGUCUUUUCUAUUAUUUCUAUCGAAAUCAU" +
                "GUCUUUACAACUAAAGUUACCUAUUGAAAUACUAAUCAAUGUGUAGCCAA" +
                "UUACUUUUUAAAGGUCUCCAUUUGAUGAAAAAACUUAAACUUGGACUUUG" +
                "AUUUUAUCAACAAGUUUUUUUACUUAAUUUUCUACGUGAAGUUUCUUAAG" +
                "UUUGAAACCGAGUUUUACUUUCUUGAAAAAAUACGCUAUACGUUUAAUUG" +
                "UCAAGACUUAAUUUUUAUUCUCGAUAACAAUUAUUAUAGCCUUUAAGAGA" +
                "ACUCCUUUAAAGAACAGAAUUUAAACUUCCAAUAUUUGAAUUAUAAAGAA" +
                "AAUUGGGUUCAAGAGAUAAUCUAGUGUAUCUCAGAAAACUUAGUUUACUU" +
                "UAUUUAAAACUAAAGGUUCCUUUAUCAUUCAUAAAAAACUAUUGGAGAUU" +
                "UUCACUUGGACUUGAAUUCGUUUAUAACCAAGGAAGUUCUAUU",
                rna["seq"]
			)

		self.assertEqual([rna["seq"].count("A"), rna["seq"].count("C"), rna["seq"].count("G"), rna["seq"].count("U")], rna["ntCount"])
		self.assertAlmostEqual(362601.870000, rna["mw"], places = 6)
		self.assertEqual("MG_001", rna["geneId"])
		self.assertEqual("MG_001_MONOMER", rna["monomerId"])
		