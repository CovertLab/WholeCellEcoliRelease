#!/usr/bin/env python

"""
Fitter

Adjust simulation parameters

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/11/2013
"""

import numpy
from wholecell.utils.constants import Constants

def fitSimulation(sim, kb):
	mc = sim.states["MoleculeCounts"]
	# tc = sim.processes["Transcription"]
	# tl = sim.processes["Translation"]
	tc_elngRate = 50 # TOKB
	tc_cellCycleLength = 1 * 3600. # TOKB
	tl_elngRate = 16 # TOKB

	idx = {}

	idx["rnaExp"] = {}

	# RNA types
	rnaTypes = dict([(x["rnaId"], x["type"]) for x in kb.genes])
	mRnaIds = [
		rna['id'] for rna in kb.rnas
		if rna["unmodifiedForm"] == None and rnaTypes[rna["id"]] in ["mRNA"]
		]
	idx["rnaExp"]["mRnas"] = numpy.array([i for i, x in enumerate(kb.rnas) if x['id'] in mRnaIds])
	idx["rnaExp"]["miscRnas"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["unmodifiedForm"] == None and rnaTypes[x["id"]] in ["miscRNA"]])
	idx["rnaExp"]["rRna23Ss"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["unmodifiedForm"] == None and x["id"] in _ids["rRna23Ss"]])
	idx["rnaExp"]["rRna16Ss"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["unmodifiedForm"] == None and x["id"] in _ids["rRna16Ss"]])
	idx["rnaExp"]["rRna5Ss"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["unmodifiedForm"] == None and x["id"] in _ids["rRna5Ss"]])
	idx["rnaExp"]["tRnas"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["unmodifiedForm"] == None and x["id"] in _ids["tRnas"]])
	idx["rnaExp"]["modified"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["unmodifiedForm"] != None])
	idx["rnaExp"]["unmodified"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["unmodifiedForm"] == None])

	idx["rnaLens"] = {}
	idx["rnaLens"]["unmodified"] = idx["rnaExp"]["unmodified"]

	idx["monExp"] = {}
	idx["monExp"]["rnap_70"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["monomerId"] != None and x["id"] in ["EG10893_RNA", "EG10894_RNA", "EG10895_RNA", "EG10896_RNA"]])

	# RNA Polymerase
	idx["rnaExp"]["rnap_70"] = numpy.array([i for i, x in enumerate(kb.rnas) if x["id"] in _ids['rnap']])

	# Feist core
	idx["FeistCore"] = {}
	idx["FeistCore"]["aminoAcids"] = numpy.array([i for i, x in enumerate(_ids["FeistCore"]) if x in _ids["aminoAcids"]])
	idx["FeistCore"]["dntps"] = numpy.array([i for i, x in enumerate(_ids["FeistCore"]) if x in _ids["dntps"]])
	idx["FeistCore"]["atp"] = numpy.array([i for i, x in enumerate(_ids["FeistCore"]) if x in _ids["atp"]])
	idx["FeistCore"]["ctp"] = numpy.array([i for i, x in enumerate(_ids["FeistCore"]) if x in _ids["ctp"]])
	idx["FeistCore"]["gtp"] = numpy.array([i for i, x in enumerate(_ids["FeistCore"]) if x in _ids["gtp"]])
	idx["FeistCore"]["utp"] = numpy.array([i for i, x in enumerate(_ids["FeistCore"]) if x in _ids["utp"]])
	idx["FeistCore"]["ntps"] = numpy.array([
		idx["FeistCore"]["atp"], idx["FeistCore"]["ctp"], idx["FeistCore"]["gtp"], idx["FeistCore"]["utp"]
	]).reshape(-1)

	idx["proteinAaCounts"] = {}
	idx["proteinAaCounts"]["notSec"] = numpy.array([i for i in xrange(21) if i != 15])

	# Change
	idx["rnaExpFracs"] = dict([(x[1], x[0]) for x in enumerate(["rRna23Ss", "rRna16Ss", "rRna5Ss", "tRnas", "mRnas"])])

	massFracRNAs = numpy.array([
		0.525,	# 23S rRNA
		0.271,	# 16S rRNA
		0.017,	# 5S rRNA
		0.146,	# tRNA
		0.041,	# mRNA (include miscRNAs here if applicable (i.e., if not setting their expression to zero))
		]) # TOKB

	mwRNAs = numpy.array([
		mc.molMass(_ids['rRna23Ss']).mean(),
		mc.molMass(_ids['rRna16Ss']).mean(),
		mc.molMass(_ids['rRna5Ss']).mean(),
		mc.molMass(_ids['tRnas']).mean(),
		# mc.molMass(_ids['matureMrnaMiscRna']).mean(), # Use if including miscRNAs in mRNA mass fraction
		mc.molMass(mRnaIds).mean()
		])

	rnaExpFracs = massFracRNAs / mwRNAs
	rnaExpFracs /= numpy.sum(rnaExpFracs)

	mc.rnaExp[idx["rnaExp"]["rRna23Ss"]] = rnaExpFracs[idx["rnaExpFracs"]["rRna23Ss"]] * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna23Ss"].size)) * numpy.ones(idx["rnaExp"]["rRna23Ss"].size)
	mc.rnaExp[idx["rnaExp"]["rRna16Ss"]] = rnaExpFracs[idx["rnaExpFracs"]["rRna16Ss"]] * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna16Ss"].size)) * numpy.ones(idx["rnaExp"]["rRna16Ss"].size)
	mc.rnaExp[idx["rnaExp"]["rRna5Ss"]]  = rnaExpFracs[idx["rnaExpFracs"]["rRna5Ss"]]  * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna5Ss"].size))  * numpy.ones(idx["rnaExp"]["rRna5Ss"].size )
	mc.rnaExp[idx["rnaExp"]["tRnas"]]    = rnaExpFracs[idx["rnaExpFracs"]["tRnas"]]    * 1. / numpy.sum(mc.rnaExp[idx["rnaExp"]["tRnas"]])          * mc.rnaExp[idx["rnaExp"]["tRnas"]]
	mc.rnaExp[idx["rnaExp"]["mRnas"]]    = rnaExpFracs[idx["rnaExpFracs"]["mRnas"]]    * 1. / numpy.sum(mc.rnaExp[idx["rnaExp"]["mRnas"]])          * mc.rnaExp[idx["rnaExp"]["mRnas"]]
	mc.rnaExp[idx["rnaExp"]["miscRnas"]] = 0. # Uncomment if producing miscRNAs
	mc.rnaExp[idx["rnaExp"]["modified"]] = 0.

	assert(numpy.abs(numpy.sum(mc.rnaExp) - 1.) < 1e-9)

	# Adjust RNAP expression levels to be the mean of the subunits
	# TODO: Account for stoichiometry
	mc.rnaExp[idx["rnaExp"]["rnap_70"]] = numpy.mean(mc.rnaExp[idx["rnaExp"]["rnap_70"]])
	mc.rnaExp /= numpy.sum(mc.rnaExp)
	mc.monExp[idx["monExp"]["rnap_70"]] = mc.rnaExp[idx["rnaExp"]["rnap_70"]] / rnaExpFracs[idx["rnaExpFracs"]["mRnas"]]
	mc.monExp /= numpy.sum(mc.monExp)

	feistCoreView = mc.countsBulkViewNew([
		(id_ + '[c]' if id_[-1] != ']' else id_) for id_ in _ids['FeistCore'] # hack to include compartments
		])

	h2oMol = mc.molecule('H2O[c]')
	h2oMass = h2oMol.massSingle()

	ppiMass = mc.molecule('PPI[c]').massSingle()

	ntpView = mc.countsBulkViewNew([
		id_ + '[c]' for id_ in _ids['ntps']
		])

	aaView = mc.countsBulkViewNew([
		id_ + '[c]' for id_ in _ids['aminoAcids']
		])

	hL = numpy.array([x["halfLife"] for x in kb.rnas if x["unmodifiedForm"] == None])
	mw_c_aas = mc.molMass(_ids['aminoAcids']) - h2oMass
	mw_c_ntps = mc.molMass(_ids['ntps']) - ppiMass
	mw_c_dntps = mc.molMass(_ids['dntps']) - ppiMass

	for iteration in xrange(5):
		# Estimate number of RNA Polymerases needed initially
		
		feistCoreView.countsBulkIs(
			numpy.round(mc.feistCoreVals * 1e-3 * Constants.nAvogadro * mc.initialDryMass) # this calculation seems wrong...
			)

		h2oMol.countBulkIs(
			(6.7e-13 / 1.36 + sim.randStream.normal(0, 15e-15)) / h2oMass * Constants.nAvogadro # not entirely sure what this does
			) # TOKB

		ntpsToPolym = numpy.round((1 - mc.fracInitFreeNTPs) * ntpView.countsBulk().sum()) # number of NTPs as RNA
		numRnas = numpy.round(ntpsToPolym / (numpy.dot(mc.rnaExp, mc.rnaLens))) # expected number of RNAs?
		
		numRnapsNeeded = numpy.sum(
			mc.rnaLens[idx["rnaLens"]["unmodified"]].astype("float") / tc_elngRate * (
				numpy.log(2) / tc_cellCycleLength + numpy.log(2) / hL
				) * numRnas * mc.rnaExp[idx["rnaExp"]["unmodified"]]
			)

		#print "numRnapsNeeded: %0.1f" % numRnapsNeeded
		
		# Estimate total number of monomers
		aasToPolym = numpy.round((1 - mc.fracInitFreeAAs) * aaView.countsBulk().sum()) # number of AAs as protein
		numMons = numpy.round(aasToPolym / (numpy.dot(mc.monExp, mc.monLens))) # expected number of proteins?

		fudge = 10000
		if numpy.min(numMons * mc.monExp[idx["monExp"]["rnap_70"]] * numpy.array([1./2, 1., 1., 1.])) < fudge * numRnapsNeeded:
			# Adjust monomer expression if necessary
			# mc.monExp[idx["monExp"]["rnap_70"]] = numpy.maximum(mc.monExp[idx["monExp"]["rnap_70"]], fudge * float(numRnapsNeeded) / numMons)
			mc.monExp[idx["monExp"]["rnap_70"]] = numpy.minimum(mc.monExp[idx["monExp"]["rnap_70"]], 2*float(numRnapsNeeded) / numMons)
			mc.monExp /= numpy.sum(mc.monExp)
			# Make corresponding change to mc.rnaExp
			mc.rnaExp[idx["rnaExp"]["rnap_70"]] = rnaExpFracs[idx["rnaExpFracs"]["mRnas"]] * mc.monExp[idx["monExp"]["rnap_70"]]
			# mc.rnaExp[idx["rnaExp"]["mRnas"]] = rnaExpFracs[idx["rnaExpFracs"]["mRnas"]] * mc.rnaExp[idx["rnaExp"]["mRnas"]] / numpy.sum(mc.rnaExp[idx["rnaExp"]["mRnas"]])
			mc.rnaExp /= numpy.sum(mc.rnaExp)

		# Estimate number of ribosomes needed initially
		#numRibsNeeded = numpy.sum(mc.monLens.astype("float") / tl.elngRate * ( numpy.log(2) / tc_cellCycleLength) * numMons * mc.monExp)
		numRibsNeeded = numpy.sum(mc.monLens.astype("float") / tl_elngRate * ( numpy.log(2) / tc_cellCycleLength) * numMons * mc.monExp)
		#print "numRibsNeeded: %0.1f" % numRibsNeeded
		fudge = 1.1
		if numpy.sum(numRnas * mc.rnaExp[idx["rnaExp"]["rRna23Ss"]]) < fudge * numRibsNeeded:
			mc.rnaExp[idx["rnaExp"]["rRna23Ss"]] = numpy.maximum(mc.rnaExp[idx["rnaExp"]["rRna23Ss"]], fudge * float(numRibsNeeded) / numRnas)
			mc.rnaExp /= numpy.sum(mc.rnaExp)
			raise Exception, "Changing RNA mass fractions. Write code to handle this."
		if numpy.sum(numRnas * mc.rnaExp[idx["rnaExp"]["rRna16Ss"]]) < fudge * numRibsNeeded:
			mc.rnaExp[idx["rnaExp"]["rRna16Ss"]] = numpy.maximum(mc.rnaExp[idx["rnaExp"]["rRna16Ss"]], fudge * float(numRibsNeeded) / numRnas)
			mc.rnaExp /= numpy.sum(mc.rnaExp)
			raise Exception, "Changing RNA mass fractions. Write code to handle this."
		if numpy.sum(numRnas * mc.rnaExp[idx["rnaExp"]["rRna16Ss"]]) < fudge * numRibsNeeded:
			mc.rnaExp[idx["rnaExp"]["rRna16Ss"]] = numpy.maximum(mc.rnaExp[idx["rnaExp"]["rRna16Ss"]], fudge * float(numRibsNeeded) / numRnas)
			mc.rnaExp /= numpy.sum(mc.rnaExp)
			raise Exception, "Changing RNA mass fractions. Write code to handle this."

		# Calculate RNA Synthesis probabilities
		if 'Transcription' in sim.processes:
			tc = sim.processes['Transcription']

			hLfull = numpy.array([x["halfLife"] if x["unmodifiedForm"] == None else numpy.inf for x in kb.rnas])
			# tc.rnaSynthProb = mc.rnaLens.astype("float") / tc_elngRate * ( numpy.log(2) / tc_cellCycleLength + numpy.log(2) / hLfull ) * numRnas * mc.rnaExp
			tc.rnaSynthProb = ( numpy.log(2) / tc_cellCycleLength + numpy.log(2) / hLfull ) * numRnas * mc.rnaExp
			tc.rnaSynthProb /= numpy.sum(tc.rnaSynthProb)

		# Assert relationship between mc.monExp and mc.rnaExp
		assert(numpy.all((mc.rnaExp[idx["rnaExp"]["mRnas"]] - rnaExpFracs[idx["rnaExpFracs"]["mRnas"]] * mc.monExp) < 1e-5))

		# Align biomass with process usages
		valsOrig = mc.feistCoreVals.copy()

		normalize = lambda x: numpy.array(x).astype("float") / numpy.linalg.norm(x, 1)

		# Amino acids (Protein)
		#f_w = normalize(numpy.sum(mc.monExp.reshape(-1, 1) * tl.proteinAaCounts[:, idx["proteinAaCounts"]["notSec"]], axis = 0))
		f_w = numpy.array([ 0.09832716,  0.05611487,  0.04021716,  0.0545386 ,  0.00908125,
    						0.06433478,  0.04242188,  0.07794587,  0.02055925,  0.05964359,
					        0.09432389,  0.05520678,  0.02730249,  0.03564025,  0.04069936,
					        0.05387673,  0.05485896,  0.01133458,  0.02679389,  0.07677868]) # TOKB
		mc.feistCoreVals[idx["FeistCore"]["aminoAcids"]] = 1000 * 0.5794 * f_w / mw_c_aas # TOKB

		# NTPs (RNA)
		# f_w = numpy.array([ 0.25375551,  0.23228423,  0.30245459,  0.21150567])
		f_w = numpy.array([ 0.248,  0.238,  0.300,  0.214 ]) # TOKB
		# f_w = normalize(numpy.sum(tc.rnaSynthProb.reshape(-1, 1) * tc.rnaNtCounts, axis = 0))
		mc.feistCoreVals[idx["FeistCore"]["ntps"]] = 1000 * 0.216 * f_w / mw_c_ntps # TOKB

		# dNTPS (DNA)
		f_w = normalize(numpy.array([kb.genomeSeq.count("A"), kb.genomeSeq.count("C"), kb.genomeSeq.count("G"), kb.genomeSeq.count("T")]))
		mc.feistCoreVals[idx["FeistCore"]["dntps"]] = 1000 * 0.0327 * f_w / mw_c_dntps # TOKB

		#print "||delta biomass||_1: %0.3f" % numpy.linalg.norm(valsOrig - mc.vals["FeistCore"], 1)

	# import ipdb
	# ipdb.set_trace()

	sim.calcInitialConditions() # Recalculate initial conditions based on fit parameters

_ids = {} # TOKB
_ids["tRnas"] = [
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

_ids["rRna23Ss"] = [
	"RRLA-RRNA", "RRLB-RRNA", "RRLC-RRNA", "RRLD-RRNA", "RRLE-RRNA", "RRLG-RRNA", "RRLH-RRNA",
	]

_ids["rRna16Ss"] = [
	"RRSA-RRNA", "RRSB-RRNA", "RRSC-RRNA", "RRSD-RRNA", "RRSE-RRNA", "RRSG-RRNA", "RRSH-RRNA",
	]

_ids["rRna5Ss"]  = [
	"RRFA-RRNA", "RRFB-RRNA", "RRFC-RRNA", "RRFD-RRNA", "RRFE-RRNA", "RRFF-RRNA", "RRFG-RRNA", "RRFH-RRNA"
	]

_ids["FeistCore"] = [
	"ALA-L", "ARG-L", "ASN-L", "ASP-L", "CYS-L", "GLN-L", "GLU-L", "GLY", "HIS-L", "ILE-L",
	"LEU-L", "LYS-L", "MET-L", "PHE-L", "PRO-L", "SER-L", "THR-L", "TRP-L", "TYR-L", "VAL-L",
	"DATP", "DCTP", "DGTP", "DTTP", "CTP", "GTP", "UTP", "ATP", "MUREIN5PX4P[p]", "KDO2LIPID4[o]",
	"PE160", "PE161", "K", "NH4", "MG2", "CA2", "FE2", "FE3", "CU2", "MN2",
	"MOBD", "COBALT2", "ZN2", "CL", "SO4", "PI", "COA", "NAD", "NADP", "FAD",
	"THF", "MLTHF", "10FTHF", "THMPP", "PYDX5P", "PHEME", "SHEME", "UDCPDP", "AMET", "2OHPH",
	"RIBFLV"
	]

_ids["aminoAcids"] = [
	"ALA-L", "ARG-L", "ASN-L", "ASP-L", "CYS-L", "GLN-L", "GLU-L", "GLY", "HIS-L", "ILE-L",
	"LEU-L", "LYS-L", "MET-L", "PHE-L", "PRO-L", "SER-L", "THR-L", "TRP-L", "TYR-L", "VAL-L",
	]

_ids["ntps"] = [
	"ATP", "CTP", "GTP", "UTP"
	]

_ids["dntps"] = [
	"DATP", "DCTP", "DGTP", "DTTP"
	]

_ids["atp"], _ids["ctp"], _ids["gtp"], _ids["utp"] = (
	["ATP"], ["CTP"], ["GTP"], ["UTP"]
	)

_ids['rnap'] = ["EG10893_RNA", "EG10894_RNA", "EG10895_RNA", "EG10896_RNA"]
