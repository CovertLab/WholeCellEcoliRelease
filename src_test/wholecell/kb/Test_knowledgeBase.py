"""
Test KnowledgeBase.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/22/2013
"""

import unittest
import warnings

import numpy
import cPickle
import os
import wholecell.kb.KnowledgeBase

class Test_randStream(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.kb = cPickle.load(open(os.path.join("data", "fixtures", "KnowledgeBase.cPickle"), "r"))
		
		# To load from the "raw" data, uncomment the following:
		# self.kb = wholecell.kb.KnowledgeBase.KnowledgeBase(dataFileName = "data/KnowledgeBase.xlsx",
		# 													 seqFileName = "data/KnowledgeBase.fna")

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
		self.assertNotEqual(None, met)
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
		self.assertNotEqual(None, gene)
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
		self.assertNotEqual(None, rna)
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
		self.assertTrue(numpy.array_equal([rna["seq"].count("A"), rna["seq"].count("C"), rna["seq"].count("G"), rna["seq"].count("U")], rna["ntCount"]))
		self.assertAlmostEqual(362601.870000, rna["mw"], places = 6)
		self.assertEqual("MG_001", rna["geneId"])
		self.assertEqual("MG_001_MONOMER", rna["monomerId"])

	def test_proteins(self):
		kb = self.kb

		self.assertEqual(482 + 164, len(kb.proteins))
		self.assertEqual(482, len([x for x in kb.proteins if x["monomer"] == True]))
		self.assertEqual(164, len([x for x in kb.proteins if x["monomer"] == False]))

		# Monomers
		mon = next((x for x in kb.proteins if x["id"] == "MG_001_MONOMER"), None)
		self.assertNotEqual(None, mon)
		self.assertTrue(dict, type(mon))
		self.assertEqual("DNA polymerase III, beta subunit", mon["name"])
		self.assertTrue(mon["monomer"])
		self.assertEqual(0, len(mon["composition"]))
		self.assertEqual("c", mon["compartment"])
		self.assertEqual(0, len(mon["formationProcess"]))
		self.assertEqual(
	            "MKILINKSELNKILKKMNNVIISNNKIKPHHSYFLIEAKEKEINFYANNE" +
                "YFSVKCNLNKNIDILEQGSLIVKGKIFNDLINGIKEEIITIQEKDQTLLV" +
                "KTKKTSINLNTINVNEFPRIRFNEKNDLSEFNQFKINYSLLVKGIKKIFH" +
                "SVSNNREISSKFNGVNFNGSNGKEIFLEASDTYKLSVFEIKQETEPFDFI" +
                "LESNLLSFINSFNPEEDKSIVFYYRKDNKDSFSTEMLISMDNFMISYTSV" +
                "NEKFPEVNYFFEFEPETKIVVQKNELKDALQRIQTLAQNERTFLCDMQIN" +
                "SSELKIRAIVNNIGNSLEEISCLKFEGYKLNISFNPSSLLDHIESFESNE" +
                "INFDFQGNSKYFLITSKSEPELKQILVPSR",
                mon["seq"]
			)
		self.assertTrue(
			numpy.array_equal(
				[mon["seq"].count("A"), mon["seq"].count("R"), mon["seq"].count("N"), mon["seq"].count("D"), mon["seq"].count("C"),
				mon["seq"].count("E"), mon["seq"].count("Q"), mon["seq"].count("G"), mon["seq"].count("H"), mon["seq"].count("I"),
				mon["seq"].count("L"), mon["seq"].count("K"), mon["seq"].count("M"), mon["seq"].count("F"), mon["seq"].count("P"),
				mon["seq"].count("S"), mon["seq"].count("T"), mon["seq"].count("W"), mon["seq"].count("Y"), mon["seq"].count("V")],
				mon["aaCount"])
			)
		self.assertAlmostEqual(44317.8348 - (177.22 - 149.21), mon["mw"], delta = 1e-4 * mon["mw"])
		self.assertEqual("MG_001", mon["geneId"])
		self.assertEqual("MG_001", mon["rnaId"])

		# Complexes
		cpx = next((x for x in kb.proteins if x["id"] == "DNA_GYRASE"), None)
		self.assertNotEqual(None, cpx)
		self.assertTrue(dict, type(cpx))
		self.assertEqual("DNA_GYRASE", cpx["id"])
		self.assertEqual("DNA gyrase", cpx["name"])
		self.assertFalse(cpx["monomer"])
		self.assertEqual([
			{"molecule": "MG_003_MONOMER", "compartment": "c", "coeff": -2, "form": "mature", "type": "protein"},
			{"molecule": "MG_004_MONOMER", "compartment": "c", "coeff": -2, "form": "mature", "type": "protein"},
			{"molecule": "DNA_GYRASE", "compartment": "c", "coeff": 1, "form": "mature", "type": "protein"}
			], cpx["composition"]
			)
		self.assertEqual("c", cpx["compartment"])
		self.assertEqual("Complexation", cpx["formationProcess"])
		self.assertEqual("", cpx["seq"])
		self.assertTrue(
			numpy.array_equal(
				2 * next((x for x in kb.proteins if x["id"] == "MG_003_MONOMER"), None)["aaCount"] +
				2 * next((x for x in kb.proteins if x["id"] == "MG_004_MONOMER"), None)["aaCount"],
				cpx["aaCount"])
			)

		self.assertTrue(numpy.array_equal(numpy.zeros(4), cpx["ntCount"]))
		self.assertAlmostEqual(334028.2216, cpx["mw"], delta = 1e-2 * cpx["mw"])
		self.assertEqual("", cpx["geneId"])
		self.assertEqual("", cpx["rnaId"])

		cpx = next((x for x in kb.proteins if x["id"] == "MG_014_015_DIMER"), None)
		self.assertNotEqual(None, cpx)
		self.assertEqual("m", cpx["compartment"])

	def test_reactions(self):
		kb = self.kb

		self.assertEqual(643, len(kb.reactions))

		rxn = next((x for x in kb.reactions if x["id"] == "AckA"), None)
		self.assertNotEqual(None, rxn)
		self.assertEqual(dict, type(rxn))
		self.assertEqual("AckA", rxn["id"])
		self.assertEqual("acetate kinase", rxn["name"])
		self.assertEqual("Metabolism", rxn["process"])
		self.assertEqual("2.7.2.1", rxn["ec"])
		self.assertEqual(0, rxn["dir"])
		self.assertEqual([
			{"molecule": "ACTP", "form": "mature", "compartment": "c", "coeff": -1, "type": "metabolite"},
			{"molecule": "ADP", "form": "mature", "compartment": "c", "coeff": -1, "type": "metabolite"},
			{"molecule": "AC", "form": "mature", "compartment": "c", "coeff": 1, "type": "metabolite"},
			{"molecule": "ATP", "form": "mature", "compartment": "c", "coeff": 1, "type": "metabolite"}
			], rxn["stoichiometry"]
			)
		self.assertEqual(
			{"id": "MG_357_DIMER", "form": "mature", "compartment": "c",
			 "kCatFor": 68.0 / 60.0 * 1e-3 * next((x for x in kb.proteins if x["id"] == "MG_357_DIMER"), None)["mw"],
			 "kCatRev": 70.0 / 60.0 * 1e-3 * next((x for x in kb.proteins if x["id"] == "MG_357_DIMER"), None)["mw"]
			}, rxn["enzyme"]
			)

		rxn = next((x for x in kb.reactions if x["id"] == "Aas1"), None)
		self.assertEqual(1, rxn["dir"])

		rxn = next((x for x in kb.reactions if x["id"] == "Cls1"), None)
		self.assertEqual([
			{"molecule": "PG160", "form": "mature", "compartment": "m", "coeff": -2, "type": "metabolite"},
			{"molecule": "CL160", "form": "mature", "compartment": "m", "coeff": 1, "type": "metabolite"},
			{"molecule": "GL", "form": "mature", "compartment": "c", "coeff": 1, "type": "metabolite"}
			], rxn["stoichiometry"]
			)