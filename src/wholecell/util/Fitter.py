#!/usr/bin/env python

"""
Fitter

Adjust simulation parameters

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/11/2013
"""

import numpy

class Fitter(object):
	""" Fitter """

	@staticmethod
	def FitSimulation(sim, kb):
		mc = sim.getState("MoleculeCounts")
		tc = sim.getProcess("Transcription")
		tl = sim.getProcess("Translation")

		ids = {}
		ids["tRnas"] = [
			"gltV-tRNA", "gltT-tRNA", "gltW-tRNA", "gltU-tRNA", "glnU-tRNA", "glnW-tRNA", "glnX-tRNA", "glnV-tRNA", "serT-tRNA", "serW-tRNA", "selC-tRNA",
			"serU-tRNA", "serV-tRNA", "serX-tRNA", "RNA0-302", "lysV-tRNA", "RNA0-303", "RNA0-301", "lysW-tRNA", "lysT-tRNA", "RNA0-306", "metY-tRNA",
			"metW-tRNA", "metZ-tRNA", "metU-tRNA", "metT-tRNA", "thrW-tRNA", "thrV-tRNA", "thrU-tRNA", "thrT-tRNA", "trpT-tRNA", "pheV-tRNA",
			"pheU-tRNA", "glyV-tRNA", "glyY-tRNA", "glyU-tRNA", "glyT-tRNA", "glyX-tRNA", "glyW-tRNA", "proL-tRNA", "proK-tRNA", "proM-tRNA",
			"RNA0-300", "valU-tRNA", "valV-tRNA", "valX-tRNA", "valY-tRNA", "valT-tRNA", "valW-tRNA", "hisR-tRNA", "ileX-tRNA", "RNA0-305",
			"ileV-tRNA", "ileT-tRNA", "ileU-tRNA", "tyrV-tRNA", "tyrU-tRNA", "tyrT-tRNA", "alaX-tRNA", "alaW-tRNA", "alaT-tRNA", "alaV-tRNA",
			"alaU-tRNA", "argY-tRNA", "argZ-tRNA", "argX-tRNA", "argU-tRNA", "argV-tRNA", "argQ-tRNA", "argW-tRNA", "aspV-tRNA", "aspU-tRNA",
			"aspT-tRNA", "RNA0-304", "asnV-tRNA", "asnU-tRNA", "asnT-tRNA", "leuU-tRNA", "leuQ-tRNA", "leuX-tRNA", "leuV-tRNA", "leuT-tRNA",
			"leuZ-tRNA", "leuW-tRNA", "leuP-tRNA", "cysT-tRNA"
		]

		ids["rRna23Ss"] = [
			"RRLA-RRNA", "RRLB-RRNA", "RRLC-RRNA", "RRLD-RRNA", "RRLE-RRNA", "RRLG-RRNA", "RRLH-RRNA",
		]
		ids["rRna16Ss"] = [
			"RRSA-RRNA", "RRSB-RRNA", "RRSC-RRNA", "RRSD-RRNA", "RRSE-RRNA", "RRSG-RRNA", "RRSH-RRNA",
		]
		ids["rRna5Ss"]  = [
			"RRFA-RRNA", "RRFB-RRNA", "RRFC-RRNA", "RRFD-RRNA", "RRFE-RRNA", "RRFF-RRNA", "RRFG-RRNA", "RRFH-RRNA"
		]

		idx = {}

		idx["rnaExp"] = {}

		# RNA types
		rnaTypes = dict([(x["rnaId"], x["type"]) for x in kb.genes])
		idx["rnaExp"]["mRnas"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["modifiedForm"] == False and rnaTypes[x["id"]] in ["mRNA", "miscRNA"]])
		idx["rnaExp"]["rRna23Ss"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["modifiedForm"] == False and x["id"] in ids["rRna23Ss"]])
		idx["rnaExp"]["rRna16Ss"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["modifiedForm"] == False and x["id"] in ids["rRna16Ss"]])
		idx["rnaExp"]["rRna5Ss"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["modifiedForm"] == False and x["id"] in ids["rRna5Ss"]])
		idx["rnaExp"]["tRnas"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["modifiedForm"] == False and x["id"] in ids["tRnas"]])
		idx["rnaExp"]["modified"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["modifiedForm"] == True])
		idx["rnaExp"]["unmodified"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["modifiedForm"] == False])

		idx["rnaLens"] = {}
		idx["rnaLens"]["unmodified"] = idx["rnaExp"]["unmodified"]

		idx["monExp"] = {}
		idx["monExp"]["rnap"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["monomerId"] != None and x["id"] in ["EG10893_RNA", "EG10894_RNA", "EG10895_RNA", "EG10896_RNA"]])

		# RNA Polymerase
		rnap_ids = ["EG10893_RNA", "EG10894_RNA", "EG10895_RNA", "EG10896_RNA"]
		idx["rnaExp"]["rnap_70"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["id"] in rnap_ids])


		# Change
		idx["rnaExpFracs"] = dict([(x[1], x[0]) for x in enumerate(["rRna23Ss", "rRna16Ss", "rRna5Ss", "tRnas", "mRnas"])])

		massFracRNAs = numpy.array([
			0.525,	# 23S rRNA
			0.271,	# 16S rRNA
			0.017,	# 5S rRNA
			0.146,	# tRNA
			0.041,	# mRNA (include miscRNAs here)
		])
		mwRNAs = numpy.array([
			numpy.mean(mc.mws[mc.idx["rRna23Ss"]]),
			numpy.mean(mc.mws[mc.idx["rRna16Ss"]]),
			numpy.mean(mc.mws[mc.idx["rRna5Ss"]]),
			numpy.mean(mc.mws[mc.idx["tRnas"]]),
			numpy.mean(mc.mws[mc.idx["matureMrnaMiscRna"]]),
		])
		rnaExpFracs = massFracRNAs / mwRNAs
		rnaExpFracs /= numpy.sum(rnaExpFracs)

		mc.rnaExp[idx["rnaExp"]["rRna23Ss"]] = rnaExpFracs[idx["rnaExpFracs"]["rRna23Ss"]] * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna23Ss"].size)) * numpy.ones(idx["rnaExp"]["rRna23Ss"].size)
		mc.rnaExp[idx["rnaExp"]["rRna16Ss"]] = rnaExpFracs[idx["rnaExpFracs"]["rRna16Ss"]] * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna16Ss"].size)) * numpy.ones(idx["rnaExp"]["rRna16Ss"].size)
		mc.rnaExp[idx["rnaExp"]["rRna5Ss"]]  = rnaExpFracs[idx["rnaExpFracs"]["rRna5Ss"]]  * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna5Ss"].size))  * numpy.ones(idx["rnaExp"]["rRna5Ss"].size )
		mc.rnaExp[idx["rnaExp"]["tRnas"]]    = rnaExpFracs[idx["rnaExpFracs"]["tRnas"]]    * 1. / numpy.sum(mc.rnaExp[idx["rnaExp"]["tRnas"]])          * mc.rnaExp[idx["rnaExp"]["tRnas"]]
		mc.rnaExp[idx["rnaExp"]["mRnas"]]    = rnaExpFracs[idx["rnaExpFracs"]["mRnas"]]    * 1. / numpy.sum(mc.rnaExp[idx["rnaExp"]["mRnas"]])          * mc.rnaExp[idx["rnaExp"]["mRnas"]]
		mc.rnaExp[idx["rnaExp"]["modified"]] = 0.

		assert(numpy.abs(numpy.sum(mc.rnaExp) - 1.) < 1e-9)

		mc.rnaExp[idx["rnaExp"]["rnap_70"]] = numpy.mean(mc.rnaExp[idx["rnaExp"]["rnap_70"]])
		mc.rnaExp /= numpy.sum(mc.rnaExp)

		# Estimate number of RNA Polymerases needed
		from wholecell.util.Constants import Constants
		mc.counts[mc.idx["FeistCoreRows"], mc.idx["FeistCoreCols"]] = numpy.round(mc.vals["FeistCore"] * 1e-3 * Constants.nAvogadro * mc.initialDryMass)
		mc.counts[mc.idx["h2o"], mc.cIdx["c"]] = (6.7e-13 / 1.36 + numpy.random.normal(0, 15e-15)) / mc.mws[mc.idx["h2o"]] * Constants.nAvogadro

		ntpsToPolym = numpy.round((1 - mc.fracInitFreeNTPs) * numpy.sum(mc.counts[mc.idx["ntps"], mc.cIdx["c"]]))
		numRnas = numpy.round(ntpsToPolym / (numpy.dot(mc.rnaExp, mc.rnaLens)))

		hL = numpy.array([x["halfLife"] for x in kb.rnas if x["modifiedForm"] == False])
		numRnapsNeeded = numpy.sum(mc.rnaLens[idx["rnaLens"]["unmodified"]].astype("float") / tc.elngRate * ( numpy.log(2) / tc.cellCycleLength + numpy.log(2) / hL ) * numRnas * mc.rnaExp[idx["rnaExp"]["unmodified"]])


		aasToPolym = numpy.round((1 - mc.fracInitFreeAAs) * numpy.sum(mc.counts[mc.idx["aas"], mc.cIdx["c"]]))
		numMons = numpy.round(aasToPolym / (numpy.dot(mc.monExp, mc.monLens)))

		fudge = 1.1
		if numpy.min(numMons * mc.monExp[idx["monExp"]["rnap"]] * numpy.array([1./2, 1., 1., 1.])) < fudge * numRnapsNeeded:
			mc.monExp[idx["monExp"]["rnap"]] = numpy.maximum(mc.monExp[idx["monExp"]["rnap"]], fudge * float(numRnapsNeeded) / numMons)
			mc.monExp /= numpy.sum(mc.monExp)

		# Estimate number of ribosomes needed
		numRibsNeeded = numpy.sum(mc.monLens.astype("float") / tl.elngRate * ( numpy.log(2) / tc.cellCycleLength) * numMons * mc.monExp)
		fudge = 1.1
		if numpy.min(numRnas * mc.rnaExp[idx["rnaExp"]["rRna23Ss"]]) < fudge * numRibsNeeded:
			mc.rnaExp[idx["rnaExp"]["rRna23Ss"]] = numpy.maximum(mc.rnaExp[idx["rnaExp"]["rRna23Ss"]], fudge * float(numRibsNeeded) / numRnas)
			mc.rnaExp /= numpy.sum(mc.rnaExp)
		if numpy.min(numRnas * mc.rnaExp[idx["rnaExp"]["rRna16Ss"]]) < fudge * numRibsNeeded:
			mc.rnaExp[idx["rnaExp"]["rRna16Ss"]] = numpy.maximum(mc.rnaExp[idx["rnaExp"]["rRna16Ss"]], fudge * float(numRibsNeeded) / numRnas)
			mc.rnaExp /= numpy.sum(mc.rnaExp)
		if numpy.min(numRnas * mc.rnaExp[idx["rnaExp"]["rRna16Ss"]]) < fudge * numRibsNeeded:
			mc.rnaExp[idx["rnaExp"]["rRna16Ss"]] = numpy.maximum(mc.rnaExp[idx["rnaExp"]["rRna16Ss"]], fudge * float(numRibsNeeded) / numRnas)
			mc.rnaExp /= numpy.sum(mc.rnaExp)

		# Calculate RNA Synthesis probabilities
		hLfull = numpy.array([x["halfLife"] if x["modifiedForm"] == False else numpy.inf for x in kb.rnas])
		tc.rnaSynthProb = mc.rnaLens.astype("float") / tc.elngRate * ( numpy.log(2) / tc.cellCycleLength + numpy.log(2) / hLfull ) * numRnas * mc.rnaExp
		tc.rnaSynthProb /= numpy.sum(tc.rnaSynthProb)


		# idx["monExp"] = {}
		# idx["monExp"]["rnap_70"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["monomerId"] != None and x["id"] in rnap_ids])
		# mc.monExp[idx["monExp"]["rnap_70"]] = numpy.mean(mc.monExp[idx["monExp"]["rnap_70"]])
		# mc.monExp /= numpy.sum(mc.monExp)