"""
Test Simulation.py
Tests whole-cell simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/8/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib

import numpy
import cPickle
import os
#import matplotlib
#matplotlib.use("agg")
import wholecell.kb.KnowledgeBase

class Test_Simulation(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
		self.kb = cPickle.load(open(os.path.join("data","fixtures","KnowledgeBase.cPickle"), "r"))

	def tearDown(self):
		pass


	# --- Tests for run-time errors ---
	@noseAttrib.attr('smalltest')
	def test_construction(self):
		import wholecell.sim.Simulation

		# Construct simulation
		sim = wholecell.sim.Simulation.Simulation()

	@noseAttrib.attr('smalltest')
	def test_initialize(self):
		import wholecell.sim.Simulation

		# Construct simulation
		sim = wholecell.sim.Simulation.Simulation()
		sim.initialize(self.kb)

	@noseAttrib.attr('smalltest')
	def test_run(self):

		# Simulate
		sim = self.sim
		sim.setOptions({"lengthSec": 10})
		sim.run()

		self.assertEqual(10, sim.states["Time"].value)

	@noseAttrib.attr('smalltest')
	def test_disk_and_shell_logger(self):
		import wholecell.sim.logger.Shell
		import wholecell.sim.logger.Disk

		# Output directory
		outDir = os.path.join("out", "test", "SimulationTest_testLogging")

		# Run simulation
		sim = self.sim
		sim.setOptions({"lengthSec": 10})
		sim.loggerAdd(wholecell.sim.logger.Shell.Shell())
		sim.loggerAdd(wholecell.sim.logger.Disk.Disk(outDir = outDir))
		sim.run()
		

		# TODO: Finish - call from Simulation.Simulation.loadSimulation
		readPath = os.path.join(outDir, 'state.hdf')
		reloadedSim = wholecell.sim.Simulation.Simulation.loadSimulation(self.kb, readPath, timePoint = 10)

		state_keys = sim.states.keys()
		# Need to check RandStream in another way
		state_keys.pop(state_keys.index('RandStream'))
		for state_id in state_keys:
			dynamics_keys = sim.states[state_id].getDynamics().keys()

			if 'growth' in dynamics_keys:
				# Growth calculated based on difference in time-steps will not match up
				dynamics_keys.pop(dynamics_keys.index('growth'))

			for d_key in dynamics_keys:
				if isinstance(sim.states[state_id].getDynamics()[d_key], numpy.ndarray):
					self.assertEqual(sim.states[state_id].getDynamics()[d_key].tolist(),
						reloadedSim.states[state_id].getDynamics()[d_key].tolist())
				else:
					self.assertEqual(sim.states[state_id].getDynamics()[d_key],
						reloadedSim.states[state_id].getDynamics()[d_key])
		# Check RandStream
		self.assertEqual(sim.states['RandStream'].getDynamics()['value'][1].tolist(),
						reloadedSim.states['RandStream'].getDynamics()['value'][1].tolist())

	@noseAttrib.attr('smalltest')
	def test_reload_at_later_timepoint(self):
		import wholecell.sim.logger.Shell
		import wholecell.sim.logger.Disk

		# Output directory
		outDir = os.path.join("out", "test", "SimulationTest_testLogging")

		# Run simulation
		sim = self.sim
		sim.setOptions({"lengthSec": 10})
		sim.loggerAdd(wholecell.sim.logger.Shell.Shell())
		sim.loggerAdd(wholecell.sim.logger.Disk.Disk(outDir = outDir))
		sim.run()

		# TODO: Finish - call from Simulation.Simulation.loadSimulation
		readPath = os.path.join(outDir, 'state.hdf')
		reloadedSim = wholecell.sim.Simulation.Simulation.loadSimulation(self.kb, readPath, timePoint = 5)
		reloadedSim.setOptions({"lengthSec": 10})
		reloadedSim.loggerAdd(wholecell.sim.logger.Shell.Shell())

		self.assertEqual(reloadedSim.initialStep, 5)
		self.assertEqual(reloadedSim.states['Time'].value, 5.)
		reloadedSim.run()

		state_keys = sim.states.keys()
		# Need to check RandStream in another way
		state_keys.pop(state_keys.index('RandStream'))
		for state_id in state_keys:
			dynamics_keys = sim.states[state_id].getDynamics().keys()

			if 'growth' in dynamics_keys:
				# Growth calculated based on difference in time-steps will not match up
				dynamics_keys.pop(dynamics_keys.index('growth'))

			for d_key in dynamics_keys:
				if isinstance(sim.states[state_id].getDynamics()[d_key], numpy.ndarray):
					self.assertEqual(sim.states[state_id].getDynamics()[d_key].tolist(),
						reloadedSim.states[state_id].getDynamics()[d_key].tolist())
				else:
					self.assertEqual(sim.states[state_id].getDynamics()[d_key],
						reloadedSim.states[state_id].getDynamics()[d_key])
		# Check RandStream
		self.assertEqual(sim.states['RandStream'].getDynamics()['value'][1].tolist(),
						reloadedSim.states['RandStream'].getDynamics()['value'][1].tolist())

		# Delete testing files
		import shutil
		shutil.rmtree(outDir)


	@noseAttrib.attr('smalltest')
	def test_loadSimulation_method(self):
		timepoint = 0
		with self.assertRaises(Exception) as context:
			readPath = 'test.hdf'
			wholecell.sim.Simulation.Simulation.loadSimulation(self.kb, readPath, timepoint)
		self.assertEqual(context.exception.message, 'State file specified does not exist!\n')

		with self.assertRaises(Exception) as context:
			readPath = 'test.file'
			open(readPath, 'a').close()
			wholecell.sim.Simulation.Simulation.loadSimulation(self.kb, readPath, timepoint)
			os.remove(readpath)
		self.assertEqual(context.exception.message, 'State file specified is not .hdf!\n')

	# --- Test ability to remove processes from simulation ---
	@noseAttrib.attr('smalltest')
	def test_removeProcesses(self):
		sim = wholecell.sim.Simulation.Simulation(processToInclude = ['Transcription'])
		sim.initialize(self.kb)
		self.assertEqual(['Transcription'], sim.processes.keys())

	# --- Test biology ---

	# @noseAttrib.attr("ic")
	# def test_initialConditions(self):
	# 	sim = self.sim
	# 	mass = sim.states["Mass"]
	# 	mc = sim.states["MoleculeCounts"]
	# 	met = sim.states["Metabolism"]

	# 	# Calculate initial conditions
	# 	sim.setOptions({"seed": 1})
	# 	sim.calcInitialConditions()

	# 	# Metabolite counts
	# 	self.assertTrue(numpy.all(numpy.isfinite(mc.counts)))
	# 	self.assertTrue(numpy.all(mc.counts >= 0))

	# 	# Mass
	# 	cellCompIdxs = numpy.array([mass.cIdx["c"], mass.cIdx["m"]])
	# 	self.assertTrue(numpy.allclose(13.1, numpy.sum(mass.cell[cellCompIdxs]), rtol = 1e-2))
	# 	self.assertTrue(numpy.allclose(13.1 * (0.7  + 0.3 * (1 - 0.0929 - 0.6197)), numpy.sum(mass.metabolite[cellCompIdxs]), rtol = 1e-2))
	# 	self.assertTrue(numpy.allclose(numpy.sum(mass.rna[cellCompIdxs]), 13.1 * 0.3 * 0.0929, rtol = 6e-1))
	# 	self.assertTrue(numpy.allclose(13.1 * 0.3 * 0.6197, numpy.sum(mass.protein[cellCompIdxs]), rtol = 1e-1))

	# 	# Growth
	# 	self.assertTrue(numpy.allclose(numpy.log(2) / (9 * 3600) * 3600 * 13.1, met.growth, rtol = 1e-1))

	# @noseAttrib.attr("tg")
	# def test_growth(self):
	# 	from wholecell.util.Constants import Constants
	# 	import matplotlib.pyplot as plt

	# 	sim = self.sim
	# 	met = sim.states["Metabolism"]
	# 	mc = sim.states["MoleculeCounts"]
	# 	mass = sim.states["Mass"]
	# 	tl = sim.processes["Translation"]
	# 	pm = sim.processes["ProteinMaturation"]
	# 	rm = sim.processes["RnaMaturation"]
	# 	cpx = sim.processes["Complexation"]

	# 	## Simulate
	# 	# Set options
	# 	sim.setOptions({"seed": 1, "lengthSec": 1000, "timeStepSec": 10})

	# 	# Calculate initial conditions
	# 	sim.calcInitialConditions()

	# 	# Record initial state
	# 	init = sim.getDynamics()
	# 	init["Mass"]["matureRna"]       = numpy.dot(mc.mws[mc.idx["matureRna"]], numpy.sum(mc.counts[mc.idx["matureRna"]], axis = 1)) / Constants.nAvogadro * 1e15
	# 	init["Mass"]["matureMonomers"]  = numpy.dot(mc.mws[mc.idx["matureMonomers"]], numpy.sum(mc.counts[mc.idx["matureMonomers"]], axis = 1)) / Constants.nAvogadro * 1e15
	# 	init["Mass"]["matureComplexes"] = numpy.dot(mc.mws[mc.idx["matureComplexes"]], numpy.sum(mc.counts[mc.idx["matureComplexes"]], axis = 1)) / Constants.nAvogadro * 1e15

	# 	# Simulate dynamics
	# 	time = sim.states["Time"]

	# 	ntps = numpy.zeros((4, sim.lengthSec / sim.timeStepSec))
	# 	ndps = numpy.zeros((4, sim.lengthSec / sim.timeStepSec))
	# 	nmps = numpy.zeros((4, sim.lengthSec / sim.timeStepSec))
	# 	aas = numpy.zeros((20, sim.lengthSec / sim.timeStepSec))
	# 	matureRna = numpy.zeros(sim.lengthSec / sim.timeStepSec)
	# 	matureMrna = numpy.zeros(sim.lengthSec / sim.timeStepSec)
	# 	matureMonomer = numpy.zeros(sim.lengthSec / sim.timeStepSec)
	# 	matureComplex = numpy.zeros(sim.lengthSec / sim.timeStepSec)
	# 	print " Time  Growth   AAs    GTP"
	# 	for iSec in xrange(sim.timeStepSec, sim.lengthSec + 1, sim.timeStepSec):
	# 		print "%5d  %.3f  %5d  %5d" % (iSec, met.growth, numpy.sum(mc.counts[mc.idx["aas"]]), mc.counts[mc.idx["ntps"][2], 0])

	# 		time.value = iSec
	# 		sim.evolveState()

	# 		ntps[:, iSec / sim.timeStepSec - 1] = mc.counts[mc.idx["ntps"], 0]
	# 		ndps[:, iSec / sim.timeStepSec - 1] = mc.counts[mc.idx["ndps"], 0]
	# 		nmps[:, iSec / sim.timeStepSec - 1] = mc.counts[mc.idx["nmps"], 0]
	# 		aas[:, iSec / sim.timeStepSec - 1] = mc.counts[mc.idx["aas"], 0]
	# 		matureRna[iSec / sim.timeStepSec - 1] = numpy.sum(mc.counts[numpy.unravel_index(rm.matureRna.mapping, mc.counts.shape)])
	# 		matureMrna[iSec / sim.timeStepSec - 1] = numpy.sum(mc.counts[numpy.unravel_index(tl.mrna.mapping, mc.counts.shape)])
	# 		matureMonomer[iSec / sim.timeStepSec - 1] = numpy.sum(mc.counts[numpy.unravel_index(pm.matureProteinMonomer.mapping, mc.counts.shape)])
	# 		matureComplex[iSec / sim.timeStepSec - 1] = numpy.sum(mc.counts[numpy.unravel_index(cpx.complex.mapping, mc.counts.shape)])

	# 	matureNoncodingRna = matureRna - matureMrna

	# 	## Plot
	# 	time = numpy.arange(sim.timeStepSec, sim.lengthSec + 1, sim.timeStepSec)

	# 	plt.figure(1)
	# 	plt.subplot(2, 2, 1); plt.plot(time, numpy.transpose(ntps)); plt.ylabel("NTPs"); plt.xlabel("Time (s)")
	# 	plt.subplot(2, 2, 2); plt.plot(time, numpy.transpose(ndps)); plt.ylabel("NDPs")
	# 	plt.subplot(2, 2, 3); plt.plot(time, numpy.transpose(nmps)); plt.ylabel("NMPs")
	# 	plt.subplot(2, 2, 4); plt.plot(time, numpy.transpose(aas)); plt.ylabel("AAs")
	# 	plt.show()

	# 	plt.figure(2)
	# 	plt.subplot(2, 2, 1); plt.plot(time, matureNoncodingRna); plt.ylabel("Mature non-coding RNA"); plt.xlabel("Time (s)")
	# 	plt.subplot(2, 2, 1); plt.plot(time, matureMrna); plt.ylabel("Mature mRNA")
	# 	plt.subplot(2, 2, 1); plt.plot(time, matureMonomer); plt.ylabel("Mature monomers")
	# 	plt.subplot(2, 2, 1); plt.plot(time, matureComplex); plt.ylabel("Mature complexes")
	# 	plt.show()

	# 	## Assert
	# 	expCumGrowth = numpy.exp(numpy.log(2) * sim.lengthSec / 30000.0)

	# 	# Growth
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Metabolism"]["growth"], met.growth, rtol = 2e-1))

	# 	# Mass
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Mass"]["cell"], mass.cell, rtol = 2e-1))
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Mass"]["cellDry"], mass.cellDry, rtol = 2e-1))
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Mass"]["rna"], mass.rna, rtol = 2e-1))
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Mass"]["protein"], mass.protein, rtol = 2e-1))

	# 	cellCompIdxs = numpy.array([mass.cIdx["c"], mass.cIdx["m"]])
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Mass"]["metabolite"][cellCompIdxs], mass.metabolite[cellCompIdxs], rtol = 5e-1))

	# 	# Physical counts
	# 	self.assertTrue(numpy.all(numpy.isfinite(mc.counts)))
	# 	self.assertTrue(numpy.all(mc.counts >= 0))

	# 	# RNA, protein matured and complexed
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Mass"]["matureRna"],
	# 		numpy.dot(mc.mws[mc.idx["matureRna"]], numpy.sum(mc.counts[mc.idx["matureRna"]], axis = 1)) / Constants.nAvogadro * 1e15,
	# 		rtol = 2e-1))
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Mass"]["matureMonomers"],
	# 		numpy.dot(mc.mws[mc.idx["matureMonomers"]], numpy.sum(mc.counts[mc.idx["matureMonomers"]], axis = 1)) / Constants.nAvogadro * 1e15,
	# 		rtol = 2e-1))
	# 	self.assertTrue(numpy.allclose(expCumGrowth * init["Mass"]["matureComplexes"],
	# 		numpy.dot(mc.mws[mc.idx["matureComplexes"]], numpy.sum(mc.counts[mc.idx["matureComplexes"]], axis = 1)) / Constants.nAvogadro * 1e15,
	# 		rtol = 2e-1))

	# 	# No mature RNA, protein in wrong compartments
	# 	matIdxs = numpy.where(numpy.logical_and(mc.types == mc.typeVals["rna"], mc.forms == mc.formVals["mature"]))[0]
	# 	matCompIdxs = numpy.ravel_multi_index((numpy.tile(matIdxs.reshape((-1, 1)), (1, 3)), numpy.tile(numpy.arange(3), (matIdxs.size, 1))), mc.counts.shape, order="F").reshape(-1)
	# 	self.assertFalse(numpy.any(mc.counts[numpy.unravel_index(numpy.setdiff1d(matCompIdxs, mc.idx["matureRna"]), mc.counts.shape, order = "F")]))

	# 	matIdxs = numpy.where(numpy.logical_and(mc.types == mc.typeVals["protein"], mc.forms == mc.formVals["mature"]))[0]
	# 	matAllCompIdxs = numpy.ravel_multi_index((numpy.tile(matIdxs.reshape((-1, 1)), (1, 3)), numpy.tile(numpy.arange(3), (matIdxs.size, 1))), mc.counts.shape, order="F").reshape(-1)
	# 	matCompIdxs = numpy.ravel_multi_index((matIdxs, mc.localizations[matIdxs].astype(int)), mc.counts.shape, order = "F").reshape(-1)
	# 	self.assertFalse(numpy.any(mc.counts[numpy.unravel_index(numpy.setdiff1d(matAllCompIdxs, matCompIdxs), mc.counts.shape, order = "F")]))


	# def test_metabolites(self):
	# 	kb = self.kb

	# 	self.assertEqual(722, len(kb.metabolites))

	# 	self.assertEqual(82, len([x for x in kb.metabolites if x["hydrophobic"] == True]))
	# 	self.assertEqual(640, len([x for x in kb.metabolites if x["hydrophobic"] == False]))

	# 	met = next((x for x in kb.metabolites if x["id"] == "ACCOA"), None)
	# 	self.assertNotEqual(None, met)
	# 	self.assertEqual(dict, type(met))
	# 	self.assertEqual("ACCOA", met["id"])
	# 	self.assertEqual("Acetyl-CoA", met["name"])
	# 	self.assertEqual("C23H34N7O17P3S1", met["formula"])
	# 	self.assertEqual("CC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N", met["smiles"])
	# 	self.assertEqual(-4, met["charge"])
	# 	self.assertEqual(805.538, met["mw"])
	# 	self.assertEqual(False, met["hydrophobic"])
	# 	self.assertEqual(0, met["mediaConc"])
	# 	self.assertAlmostEqual(3.5812e+03, met["biomassConc"], places = 1)
	# 	self.assertAlmostEqual(3.5812e+03, met["metabolismNewFlux"], places = 1)
	# 	self.assertEqual(0, met["metabolismRecyclingFlux"])

	# 	met = next((x for x in kb.metabolites if x["id"] == "AC"), None)
	# 	self.assertAlmostEqual(0.304758, met["mediaConc"], places = 6)
		
	# def test_geneticCode(self):
	# 	kb = self.kb

	# 	self.assertEqual(4, kb.translationTable)

	# def test_genome(self):
	# 	kb = self.kb

	# 	self.assertEqual(580076, len(kb.genomeSeq))
	# 	self.assertEqual("ACGT", "".join(sorted(list(set(kb.genomeSeq)))))

	# def test_genes(self):
	# 	kb = self.kb

	# 	self.assertEqual(525, len(kb.genes))
	# 	self.assertEqual(482, len([x for x in kb.genes if x["type"] == "mRNA"]))
	# 	self.assertEqual(3, len([x for x in kb.genes if x["type"] == "rRNA"]))
	# 	self.assertEqual(4, len([x for x in kb.genes if x["type"] == "sRNA"]))
	# 	self.assertEqual(36, len([x for x in kb.genes if x["type"] == "tRNA"]))

	# 	gene = next((x for x in kb.genes if x["id"] == "MG_001"), None)
	# 	self.assertNotEqual(None, gene)
	# 	self.assertEqual(dict, type(gene))
	# 	self.assertEqual("MG_001", gene["id"])
	# 	self.assertEqual("DNA polymerase III, beta subunit", gene["name"])
	# 	self.assertEqual("dnaN", gene["symbol"])
	# 	self.assertEqual("mRNA", gene["type"])
	# 	self.assertEqual(686, gene["start"])
	# 	self.assertEqual(1143, gene["len"])
	# 	self.assertTrue(gene["dir"])
	# 	self.assertEqual(
 #                "ATGAAAATATTAATTAATAAAAGTGAATTGAATAAAATTTTGAAAAAAAT" +
 #                "GAATAACGTTATTATTTCCAATAACAAAATAAAACCACATCATTCATATT" +
 #                "TTTTAATAGAGGCAAAAGAAAAAGAAATAAACTTTTATGCTAACAATGAA" +
 #                "TACTTTTCTGTCAAATGTAATTTAAATAAAAATATTGATATTCTTGAACA" +
 #                "AGGCTCCTTAATTGTTAAAGGAAAAATTTTTAACGATCTTATTAATGGCA" +
 #                "TAAAAGAAGAGATTATTACTATTCAAGAAAAAGATCAAACACTTTTGGTT" +
 #                "AAAACAAAAAAAACAAGTATTAATTTAAACACAATTAATGTGAATGAATT" +
 #                "TCCAAGAATAAGGTTTAATGAAAAAAACGATTTAAGTGAATTTAATCAAT" +
 #                "TCAAAATAAATTATTCACTTTTAGTAAAAGGCATTAAAAAAATTTTTCAC" +
 #                "TCAGTTTCAAATAATCGTGAAATATCTTCTAAATTTAATGGAGTAAATTT" +
 #                "CAATGGATCCAATGGAAAAGAAATATTTTTAGAAGCTTCTGACACTTATA" +
 #                "AACTATCTGTTTTTGAGATAAAGCAAGAAACAGAACCATTTGATTTCATT" +
 #                "TTGGAGAGTAATTTACTTAGTTTCATTAATTCTTTTAATCCTGAAGAAGA" +
 #                "TAAATCTATTGTTTTTTATTACAGAAAAGATAATAAAGATAGCTTTAGTA" +
 #                "CAGAAATGTTGATTTCAATGGATAACTTTATGATTAGTTACACATCGGTT" +
 #                "AATGAAAAATTTCCAGAGGTAAACTACTTTTTTGAATTTGAACCTGAAAC" +
 #                "TAAAATAGTTGTTCAAAAAAATGAATTAAAAGATGCACTTCAAAGAATTC" +
 #                "AAACTTTGGCTCAAAATGAAAGAACTTTTTTATGCGATATGCAAATTAAC" +
 #                "AGTTCTGAATTAAAAATAAGAGCTATTGTTAATAATATCGGAAATTCTCT" +
 #                "TGAGGAAATTTCTTGTCTTAAATTTGAAGGTTATAAACTTAATATTTCTT" +
 #                "TTAACCCAAGTTCTCTATTAGATCACATAGAGTCTTTTGAATCAAATGAA" +
 #                "ATAAATTTTGATTTCCAAGGAAATAGTAAGTATTTTTTGATAACCTCTAA" +
 #                "AAGTGAACCTGAACTTAAGCAAATATTGGTTCCTTCAAGATAA",
 #                gene["seq"]
	# 		)
	# 	self.assertEqual("MG_001", gene["rnaId"])

	# def test_rnas(self):
	# 	kb = self.kb

	# 	self.assertEqual(525, len(kb.rnas))

	# 	rna = next((x for x in kb.rnas if x["id"] == "MG_001"), None)
	# 	self.assertNotEqual(None, rna)
	# 	self.assertEqual(dict, type(rna))
	# 	self.assertEqual("MG_001", rna["id"])
	# 	self.assertEqual("DNA polymerase III, beta subunit", rna["name"])
	# 	self.assertEqual("mRNA", rna["type"])
	# 	self.assertAlmostEqual(8.8983e-5, rna["exp"], places = 9)
	# 	self.assertAlmostEqual(146.9388, rna["halfLife"], places = 4)
	# 	self.assertEqual(
	# 		    "UACUUUUAUAAUUAAUUAUUUUCACUUAACUUAUUUUAAAACUUUUUUUA" +
 #                "CUUAUUGCAAUAAUAAAGGUUAUUGUUUUAUUUUGGUGUAGUAAGUAUAA" +
 #                "AAAAUUAUCUCCGUUUUCUUUUUCUUUAUUUGAAAAUACGAUUGUUACUU" +
 #                "AUGAAAAGACAGUUUACAUUAAAUUUAUUUUUAUAACUAUAAGAACUUGU" +
 #                "UCCGAGGAAUUAACAAUUUCCUUUUUAAAAAUUGCUAGAAUAAUUACCGU" +
 #                "AUUUUCUUCUCUAAUAAUGAUAAGUUCUUUUUCUAGUUUGUGAAAACCAA" +
 #                "UUUUGUUUUUUUUGUUCAUAAUUAAAUUUGUGUUAAUUACACUUACUUAA" +
 #                "AGGUUCUUAUUCCAAAUUACUUUUUUUGCUAAAUUCACUUAAAUUAGUUA" +
 #                "AGUUUUAUUUAAUAAGUGAAAAUCAUUUUCCGUAAUUUUUUUAAAAAGUG" +
 #                "AGUCAAAGUUUAUUAGCACUUUAUAGAAGAUUUAAAUUACCUCAUUUAAA" +
 #                "GUUACCUAGGUUACCUUUUCUUUAUAAAAAUCUUCGAAGACUGUGAAUAU" +
 #                "UUGAUAGACAAAAACUCUAUUUCGUUCUUUGUCUUGGUAAACUAAAGUAA" +
 #                "AACCUCUCAUUAAAUGAAUCAAAGUAAUUAAGAAAAUUAGGACUUCUUCU" +
 #                "AUUUAGAUAACAAAAAAUAAUGUCUUUUCUAUUAUUUCUAUCGAAAUCAU" +
 #                "GUCUUUACAACUAAAGUUACCUAUUGAAAUACUAAUCAAUGUGUAGCCAA" +
 #                "UUACUUUUUAAAGGUCUCCAUUUGAUGAAAAAACUUAAACUUGGACUUUG" +
 #                "AUUUUAUCAACAAGUUUUUUUACUUAAUUUUCUACGUGAAGUUUCUUAAG" +
 #                "UUUGAAACCGAGUUUUACUUUCUUGAAAAAAUACGCUAUACGUUUAAUUG" +
 #                "UCAAGACUUAAUUUUUAUUCUCGAUAACAAUUAUUAUAGCCUUUAAGAGA" +
 #                "ACUCCUUUAAAGAACAGAAUUUAAACUUCCAAUAUUUGAAUUAUAAAGAA" +
 #                "AAUUGGGUUCAAGAGAUAAUCUAGUGUAUCUCAGAAAACUUAGUUUACUU" +
 #                "UAUUUAAAACUAAAGGUUCCUUUAUCAUUCAUAAAAAACUAUUGGAGAUU" +
 #                "UUCACUUGGACUUGAAUUCGUUUAUAACCAAGGAAGUUCUAUU",
 #                rna["seq"]
	# 		)
	# 	self.assertTrue(numpy.array_equal([rna["seq"].count("A"), rna["seq"].count("C"), rna["seq"].count("G"), rna["seq"].count("U")], rna["ntCount"]))
	# 	self.assertAlmostEqual(362601.870000, rna["mw"], places = 6)
	# 	self.assertEqual("MG_001", rna["geneId"])
	# 	self.assertEqual("MG_001_MONOMER", rna["monomerId"])

	# def test_proteins(self):
	# 	kb = self.kb

	# 	self.assertEqual(482 + 164, len(kb.proteins))
	# 	self.assertEqual(482, len([x for x in kb.proteins if x["monomer"] == True]))
	# 	self.assertEqual(164, len([x for x in kb.proteins if x["monomer"] == False]))

	# 	# Monomers
	# 	mon = next((x for x in kb.proteins if x["id"] == "MG_001_MONOMER"), None)
	# 	self.assertNotEqual(None, mon)
	# 	self.assertTrue(dict, type(mon))
	# 	self.assertEqual("DNA polymerase III, beta subunit", mon["name"])
	# 	self.assertTrue(mon["monomer"])
	# 	self.assertEqual(0, len(mon["composition"]))
	# 	self.assertEqual("c", mon["compartment"])
	# 	self.assertEqual(0, len(mon["formationProcess"]))
	# 	self.assertEqual(
	#             "MKILINKSELNKILKKMNNVIISNNKIKPHHSYFLIEAKEKEINFYANNE" +
 #                "YFSVKCNLNKNIDILEQGSLIVKGKIFNDLINGIKEEIITIQEKDQTLLV" +
 #                "KTKKTSINLNTINVNEFPRIRFNEKNDLSEFNQFKINYSLLVKGIKKIFH" +
 #                "SVSNNREISSKFNGVNFNGSNGKEIFLEASDTYKLSVFEIKQETEPFDFI" +
 #                "LESNLLSFINSFNPEEDKSIVFYYRKDNKDSFSTEMLISMDNFMISYTSV" +
 #                "NEKFPEVNYFFEFEPETKIVVQKNELKDALQRIQTLAQNERTFLCDMQIN" +
 #                "SSELKIRAIVNNIGNSLEEISCLKFEGYKLNISFNPSSLLDHIESFESNE" +
 #                "INFDFQGNSKYFLITSKSEPELKQILVPSR",
 #                mon["seq"]
	# 		)
	# 	self.assertTrue(
	# 		numpy.array_equal(
	# 			[mon["seq"].count("A"), mon["seq"].count("R"), mon["seq"].count("N"), mon["seq"].count("D"), mon["seq"].count("C"),
	# 			mon["seq"].count("E"), mon["seq"].count("Q"), mon["seq"].count("G"), mon["seq"].count("H"), mon["seq"].count("I"),
	# 			mon["seq"].count("L"), mon["seq"].count("K"), mon["seq"].count("M"), mon["seq"].count("F"), mon["seq"].count("P"),
	# 			mon["seq"].count("S"), mon["seq"].count("T"), mon["seq"].count("W"), mon["seq"].count("Y"), mon["seq"].count("V")],
	# 			mon["aaCount"])
	# 		)
	# 	self.assertAlmostEqual(44317.8348 - (177.22 - 149.21), mon["mw"], delta = 1e-4 * mon["mw"])
	# 	self.assertEqual("MG_001", mon["geneId"])
	# 	self.assertEqual("MG_001", mon["rnaId"])

	# 	# Complexes
	# 	cpx = next((x for x in kb.proteins if x["id"] == "DNA_GYRASE"), None)
	# 	self.assertNotEqual(None, cpx)
	# 	self.assertTrue(dict, type(cpx))
	# 	self.assertEqual("DNA_GYRASE", cpx["id"])
	# 	self.assertEqual("DNA gyrase", cpx["name"])
	# 	self.assertFalse(cpx["monomer"])
	# 	self.assertEqual([
	# 		{"molecule": "MG_003_MONOMER", "compartment": "c", "coeff": -2, "form": "mature", "type": "protein"},
	# 		{"molecule": "MG_004_MONOMER", "compartment": "c", "coeff": -2, "form": "mature", "type": "protein"},
	# 		{"molecule": "DNA_GYRASE", "compartment": "c", "coeff": 1, "form": "mature", "type": "protein"}
	# 		], cpx["composition"]
	# 		)
	# 	self.assertEqual("c", cpx["compartment"])
	# 	self.assertEqual("Complexation", cpx["formationProcess"])
	# 	self.assertEqual("", cpx["seq"])
	# 	self.assertTrue(
	# 		numpy.array_equal(
	# 			2 * next((x for x in kb.proteins if x["id"] == "MG_003_MONOMER"), None)["aaCount"] +
	# 			2 * next((x for x in kb.proteins if x["id"] == "MG_004_MONOMER"), None)["aaCount"],
	# 			cpx["aaCount"])
	# 		)

	# 	self.assertTrue(numpy.array_equal(numpy.zeros(4), cpx["ntCount"]))
	# 	self.assertAlmostEqual(334028.2216, cpx["mw"], delta = 1e-2 * cpx["mw"])
	# 	self.assertEqual("", cpx["geneId"])
	# 	self.assertEqual("", cpx["rnaId"])

	# 	cpx = next((x for x in kb.proteins if x["id"] == "MG_014_015_DIMER"), None)
	# 	self.assertNotEqual(None, cpx)
	# 	self.assertEqual("m", cpx["compartment"])

	# def test_reactions(self):
	# 	kb = self.kb

	# 	self.assertEqual(643, len(kb.reactions))

	# 	rxn = next((x for x in kb.reactions if x["id"] == "AckA"), None)
	# 	self.assertNotEqual(None, rxn)
	# 	self.assertEqual(dict, type(rxn))
	# 	self.assertEqual("AckA", rxn["id"])
	# 	self.assertEqual("acetate kinase", rxn["name"])
	# 	self.assertEqual("Metabolism", rxn["process"])
	# 	self.assertEqual("2.7.2.1", rxn["ec"])
	# 	self.assertEqual(0, rxn["dir"])
	# 	self.assertEqual([
	# 		{"molecule": "ACTP", "form": "mature", "compartment": "c", "coeff": -1, "type": "metabolite"},
	# 		{"molecule": "ADP", "form": "mature", "compartment": "c", "coeff": -1, "type": "metabolite"},
	# 		{"molecule": "AC", "form": "mature", "compartment": "c", "coeff": 1, "type": "metabolite"},
	# 		{"molecule": "ATP", "form": "mature", "compartment": "c", "coeff": 1, "type": "metabolite"}
	# 		], rxn["stoichiometry"]
	# 		)
	# 	self.assertEqual(
	# 		{"id": "MG_357_DIMER", "form": "mature", "compartment": "c",
	# 		 "kCatFor": 68.0 / 60.0 * 1e-3 * next((x for x in kb.proteins if x["id"] == "MG_357_DIMER"), None)["mw"],
	# 		 "kCatRev": 70.0 / 60.0 * 1e-3 * next((x for x in kb.proteins if x["id"] == "MG_357_DIMER"), None)["mw"]
	# 		}, rxn["enzyme"]
	# 		)

	# 	rxn = next((x for x in kb.reactions if x["id"] == "Aas1"), None)
	# 	self.assertEqual(1, rxn["dir"])

	# 	rxn = next((x for x in kb.reactions if x["id"] == "Cls1"), None)
	# 	self.assertEqual([
	# 		{"molecule": "PG160", "form": "mature", "compartment": "m", "coeff": -2, "type": "metabolite"},
	# 		{"molecule": "CL160", "form": "mature", "compartment": "m", "coeff": 1, "type": "metabolite"},
	# 		{"molecule": "GL", "form": "mature", "compartment": "c", "coeff": 1, "type": "metabolite"}
	# 		], rxn["stoichiometry"]
	# 		)
