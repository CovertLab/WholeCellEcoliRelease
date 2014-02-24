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

FEIST_CORE_VALS = numpy.array([ # TODO: This needs to go in the KB
	0.513689, 0.295792, 0.241055, 0.241055, 0.091580, 0.263160, 0.263160, 0.612638, 0.094738, 0.290529,
	0.450531, 0.343161, 0.153686, 0.185265, 0.221055, 0.215792, 0.253687, 0.056843, 0.137896, 0.423162,
	0.026166, 0.027017, 0.027017, 0.026166, 0.133508, 0.215096, 0.144104, 0.174831, 0.013894, 0.019456,
	0.063814, 0.075214, 0.177645, 0.011843, 0.007895, 0.004737, 0.007106, 0.007106, 0.003158, 0.003158,
	0.003158, 0.003158, 0.003158, 0.004737, 0.003948, 0.003948, 0.000576, 0.001831, 0.000447, 0.000223,
	0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000055, 0.000223, 0.000223,
	0.000223		# mmol/gDCW (supp info 3, "biomass_core", column G)
	]) # TOKB

INITIAL_DRY_MASS = 2.8e-13 / 1.36 # TOKB

def fitSimulation(sim, kb):
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
			numpy.array([rna['mw'] for rna in kb.rnas if rna['id'] in ids]).mean()
			for ids in (_ids['rRna23Ss'], _ids['rRna16Ss'], _ids['rRna5Ss'], _ids['tRnas'], mRnaIds)
			])

	rnaLens = numpy.array(map(lambda rna: numpy.sum(rna["ntCount"]), kb.rnas))
	rnaExpFracs = massFracRNAs / mwRNAs
	rnaExpFracs /= numpy.sum(rnaExpFracs)

	rnaExp = numpy.array([x["expression"] for x in kb.rnas])
	rnaExp /= numpy.sum(rnaExp)

	rnaExp[idx["rnaExp"]["rRna23Ss"]] = rnaExpFracs[idx["rnaExpFracs"]["rRna23Ss"]] * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna23Ss"].size)) * numpy.ones(idx["rnaExp"]["rRna23Ss"].size)
	rnaExp[idx["rnaExp"]["rRna16Ss"]] = rnaExpFracs[idx["rnaExpFracs"]["rRna16Ss"]] * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna16Ss"].size)) * numpy.ones(idx["rnaExp"]["rRna16Ss"].size)
	rnaExp[idx["rnaExp"]["rRna5Ss"]]  = rnaExpFracs[idx["rnaExpFracs"]["rRna5Ss"]]  * 1. / numpy.sum(numpy.ones(idx["rnaExp"]["rRna5Ss"].size))  * numpy.ones(idx["rnaExp"]["rRna5Ss"].size )
	rnaExp[idx["rnaExp"]["tRnas"]]    = rnaExpFracs[idx["rnaExpFracs"]["tRnas"]]    * 1. / numpy.sum(rnaExp[idx["rnaExp"]["tRnas"]])          * rnaExp[idx["rnaExp"]["tRnas"]]
	rnaExp[idx["rnaExp"]["mRnas"]]    = rnaExpFracs[idx["rnaExpFracs"]["mRnas"]]    * 1. / numpy.sum(rnaExp[idx["rnaExp"]["mRnas"]])          * rnaExp[idx["rnaExp"]["mRnas"]]
	rnaExp[idx["rnaExp"]["miscRnas"]] = 0. # Uncomment if producing miscRNAs
	rnaExp[idx["rnaExp"]["modified"]] = 0.

	assert(numpy.abs(numpy.sum(rnaExp) - 1.) < 1e-9)

	# Adjust RNAP expression levels to be the mean of the subunits
	# TODO: Account for stoichiometry
	rnaExp[idx["rnaExp"]["rnap_70"]] = numpy.mean(rnaExp[idx["rnaExp"]["rnap_70"]])
	rnaExp /= numpy.sum(rnaExp)

	mons = [x for x in kb.proteins if len(x["composition"]) == 0 and x["unmodifiedForm"] == None]
	monLens = numpy.array(map(lambda mon: numpy.sum(mon["aaCount"]), mons)) # TODO: remove map/lambda
	rnaIdToExp = dict([(x["id"], x["expression"]) for x in kb.rnas if x["monomerId"] != None])
	monExp = numpy.array([rnaIdToExp[x["rnaId"]] for x in mons])
	monExp /= numpy.sum(monExp)

	h2oMass = [met['mw7.2'] for met in kb.metabolites if met['id'] == 'H2O'][0] # is there a better way to do this?

	ppiMass = [met['mw7.2'] for met in kb.metabolites if met['id'] == 'PPI'][0]

	# TODO: separate count arrays for ntps/aas, created from feist core values

	halflife = numpy.array([x["halfLife"] for x in kb.rnas if x["unmodifiedForm"] == None])
	halflifeFull = numpy.array([x["halfLife"] if x["unmodifiedForm"] == None else numpy.inf for x in kb.rnas])
	mw_c_aas = numpy.array([met['mw7.2'] for met in kb.metabolites if met['id'] in _ids['aminoAcids']]) - h2oMass
	mw_c_ntps = numpy.array([met['mw7.2'] for met in kb.metabolites if met['id'] in _ids['ntps']]) - ppiMass
	mw_c_dntps = numpy.array([met['mw7.2'] for met in kb.metabolites if met['id'] in _ids['dntps']]) - ppiMass

	fracInitFreeNTPs = 0.0015 # TOKB
	fracInitFreeAAs = 0.001 # TOKB

	feistCoreVals = FEIST_CORE_VALS # TOKB
	initialDryMass = INITIAL_DRY_MASS # TOKB

	feistCoreCounts = numpy.round(
		feistCoreVals * 1e-3 * Constants.nAvogadro * initialDryMass
		)

	totalNTPs = feistCoreCounts[[i for i, id_ in enumerate(_ids['FeistCore']) if id_ in _ids['ntps']]].sum()
	totalAAs = feistCoreCounts[[i for i, id_ in enumerate(_ids['FeistCore']) if id_ in _ids['aminoAcids']]].sum()

	for iteration in xrange(5):
		# Estimate number of RNA Polymerases needed initially

		ntpsToPolym = numpy.round((1 - fracInitFreeNTPs) * totalNTPs) # number of NTPs as RNA
		numRnas = numpy.round(ntpsToPolym / (numpy.dot(rnaExp, rnaLens))) # expected number of RNAs?
		
		numRnapsNeeded = numpy.sum(
			rnaLens[idx["rnaLens"]["unmodified"]].astype("float") / tc_elngRate * (
				numpy.log(2) / tc_cellCycleLength + numpy.log(2) / halflife
				) * numRnas * rnaExp[idx["rnaExp"]["unmodified"]]
			)

		#print "numRnapsNeeded: %0.1f" % numRnapsNeeded
		
		# Estimate total number of monomers
		aasToPolym = numpy.round((1 - fracInitFreeAAs) * totalAAs) # number of AAs as protein
		numMons = numpy.round(aasToPolym / (numpy.dot(monExp, monLens))) # expected number of proteins?

		fudge = 10000
		if numpy.min(numMons * monExp[idx["monExp"]["rnap_70"]] * numpy.array([1./2, 1., 1., 1.])) < fudge * numRnapsNeeded:
			# Adjust monomer expression if necessary
			# monExp[idx["monExp"]["rnap_70"]] = numpy.maximum(monExp[idx["monExp"]["rnap_70"]], fudge * float(numRnapsNeeded) / numMons)
			monExp[idx["monExp"]["rnap_70"]] = numpy.minimum(monExp[idx["monExp"]["rnap_70"]], 2*float(numRnapsNeeded) / numMons)
			monExp /= numpy.sum(monExp)
			# Make corresponding change to rnaExp
			rnaExp[idx["rnaExp"]["rnap_70"]] = rnaExpFracs[idx["rnaExpFracs"]["mRnas"]] * monExp[idx["monExp"]["rnap_70"]]
			# rnaExp[idx["rnaExp"]["mRnas"]] = rnaExpFracs[idx["rnaExpFracs"]["mRnas"]] * rnaExp[idx["rnaExp"]["mRnas"]] / numpy.sum(rnaExp[idx["rnaExp"]["mRnas"]])
			rnaExp /= numpy.sum(rnaExp)

		# Estimate number of ribosomes needed initially
		#numRibsNeeded = numpy.sum(monLens.astype("float") / tl.elngRate * ( numpy.log(2) / tc_cellCycleLength) * numMons * monExp)
		numRibsNeeded = numpy.sum(monLens.astype("float") / tl_elngRate * ( numpy.log(2) / tc_cellCycleLength) * numMons * monExp)
		#print "numRibsNeeded: %0.1f" % numRibsNeeded
		fudge = 1.1
		if numpy.sum(numRnas * rnaExp[idx["rnaExp"]["rRna23Ss"]]) < fudge * numRibsNeeded:
			rnaExp[idx["rnaExp"]["rRna23Ss"]] = numpy.maximum(rnaExp[idx["rnaExp"]["rRna23Ss"]], fudge * float(numRibsNeeded) / numRnas)
			rnaExp /= numpy.sum(rnaExp)
			raise Exception, "Changing RNA mass fractions. Write code to handle this."
		if numpy.sum(numRnas * rnaExp[idx["rnaExp"]["rRna16Ss"]]) < fudge * numRibsNeeded:
			rnaExp[idx["rnaExp"]["rRna16Ss"]] = numpy.maximum(rnaExp[idx["rnaExp"]["rRna16Ss"]], fudge * float(numRibsNeeded) / numRnas)
			rnaExp /= numpy.sum(rnaExp)
			raise Exception, "Changing RNA mass fractions. Write code to handle this."
		if numpy.sum(numRnas * rnaExp[idx["rnaExp"]["rRna16Ss"]]) < fudge * numRibsNeeded:
			rnaExp[idx["rnaExp"]["rRna16Ss"]] = numpy.maximum(rnaExp[idx["rnaExp"]["rRna16Ss"]], fudge * float(numRibsNeeded) / numRnas)
			rnaExp /= numpy.sum(rnaExp)
			raise Exception, "Changing RNA mass fractions. Write code to handle this."

		# Assert relationship between monExp and rnaExp
		assert(numpy.all((rnaExp[idx["rnaExp"]["mRnas"]] - rnaExpFracs[idx["rnaExpFracs"]["mRnas"]] * monExp) < 1e-5))

	# Align biomass with process usages

	# Amino acids (Protein)
	#f_w = normalize(numpy.sum(monExp.reshape(-1, 1) * tl.proteinAaCounts[:, idx["proteinAaCounts"]["notSec"]], axis = 0))
	f_w = numpy.array([ 0.09832716,  0.05611487,  0.04021716,  0.0545386 ,  0.00908125,
						0.06433478,  0.04242188,  0.07794587,  0.02055925,  0.05964359,
						0.09432389,  0.05520678,  0.02730249,  0.03564025,  0.04069936,
						0.05387673,  0.05485896,  0.01133458,  0.02679389,  0.07677868]) # TOKB
	feistCoreVals[idx["FeistCore"]["aminoAcids"]] = 1000 * 0.5794 * f_w / mw_c_aas # TOKB

	# NTPs (RNA)
	# f_w = numpy.array([ 0.25375551,  0.23228423,  0.30245459,  0.21150567])
	f_w = numpy.array([ 0.248,  0.238,  0.300,  0.214 ]) # TOKB
	# f_w = normalize(numpy.sum(tc.rnaSynthProb.reshape(-1, 1) * tc.rnaNtCounts, axis = 0))
	feistCoreVals[idx["FeistCore"]["ntps"]] = 1000 * 0.216 * f_w / mw_c_ntps # TOKB

	# dNTPS (DNA)
	f_w = normalize(numpy.array([kb.genomeSeq.count("A"), kb.genomeSeq.count("C"), kb.genomeSeq.count("G"), kb.genomeSeq.count("T")]))
	feistCoreVals[idx["FeistCore"]["dntps"]] = 1000 * 0.0327 * f_w / mw_c_dntps # TOKB

	# NOTE: reactivate this line once we start fitting the RNA expression
	# for i, rna in enumerate(kb.rnas):
	# 	rna['expression'] = rnaExp[i]

	rnaSynthProb = ( numpy.log(2) / tc_cellCycleLength + numpy.log(2) / halflifeFull ) * numRnas * rnaExp
	rnaSynthProb /= rnaSynthProb.sum()

	# do the same for monomer expression

	sim.states['MoleculeCounts'].rnaExp = rnaExp
	sim.states['MoleculeCounts'].monExp = monExp
	sim.states['MoleculeCounts'].feistCoreVals = feistCoreVals
	
	# Calculate RNA Synthesis probabilities
	if 'Transcription' in sim.processes:
		sim.processes['Transcription'].rnaSynthProb = rnaSynthProb

	sim.calcInitialConditions() # Recalculate initial conditions based on fit parameters

	# TODO: return/save fitted KB instead of a modified simulation

def normalize(array):
	return numpy.array(array).astype("float") / numpy.linalg.norm(array, 1)

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
