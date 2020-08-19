from __future__ import absolute_import, division, print_function

import os
import yaml
import json
import numpy as np

from reconstruction.spreadsheets import read_tsv


METABOLITE_MASS_FILE = os.path.join("reconstruction", "ecoli", "flat", "metabolites.tsv")
MONOMER_MASS_FILE = os.path.join("reconstruction", "ecoli", "flat", "proteins.tsv")
RNA_MASS_FILE = os.path.join("reconstruction", "ecoli", "flat", "rnas.tsv")
WATER_MASS_FILE = os.path.join("reconstruction", "ecoli", "flat", "water.tsv")
EXISTING_PROTEIN_COMPLEX_FILE = os.path.join("reconstruction", "ecoli", "flat", "proteinComplexes.tsv")
NEW_PROTEIN_COMPLEX_FILE = os.path.join("reconstruction", "ecoli", "flat", "proteinComplexes_new.tsv")
NEW_COMPLEXATION_REACTION_FILE = os.path.join("reconstruction", "ecoli", "flat", "complexationReactions_new.tsv")
NEW_EQUILIBRIUM_REACTION_FILE = os.path.join("reconstruction", "ecoli", "flat", "equilibriumReactions_new.tsv")
EXISTING_EQUILIBRIUM_REACTION_FILE = os.path.join("reconstruction", "ecoli", "flat", "equilibriumReactions.tsv")

ECOCYC_DUMP = os.path.join("reconstruction", "ecoli", "flat", "eco_wc_test_fun.json")

DEFAULT_EQUILIBRIUM_BINDING_RATES = (1, 1e-06)

REACTION_ID_BLACKLIST = [
	"CPLX0-3964",	# Full ribosome--this gets formed in the simulation
	"RNAP32-CPLX",	# Full RNA Polymerase--simulation only uses apo form
	"RNAP54-CPLX",	# Full RNA Polymerase--simulation only uses apo form
	"RNAP70-CPLX",	# Full RNA Polymerase--simulation only uses apo form
	"RNAPE-CPLX",	# Full RNA Polymerase--simulation only uses apo form
	"RNAPS-CPLX",	# Full RNA Polymerase--simulation only uses apo form
	"CPLX0-221",	# Full RNA Polymerase--simulation only uses apo form
	"CPLX0-222",	# Full RNA Polymerase--simulation only uses apo form
]

def getMetaboliteMasses():
	data = read_tsv(METABOLITE_MASS_FILE)
	D = dict([(x["id"], np.array([0, 0, 0, 0, 0, 0, 0, x["mw"], 0, 0, 0])) for x in data])
	data = read_tsv(WATER_MASS_FILE)
	E = dict([(x["id"], np.array([0, 0, 0, 0, 0, 0, 0, 0, x["mw"], 0, 0])) for x in data])
	return dict(D, **E)

def getMonomerMasses():
	data = read_tsv(MONOMER_MASS_FILE)
	D = dict([(x["id"], np.array(x["mw"])) for x in data])
	return D

def getRnaMasses():
	data = read_tsv(RNA_MASS_FILE)
	D = dict([(x["id"], np.array(x["mw"])) for x in data])
	return D

def getNames():
	data = read_tsv(EXISTING_PROTEIN_COMPLEX_FILE)
	D = dict([(x["id"], x["name"]) for x in data])
	return D

def getEquilibriumBindingRates():
	data = read_tsv(EXISTING_EQUILIBRIUM_REACTION_FILE)
	D = dict((x["id"], (x["forward rate"], x["reverse rate"])) for x in data)
	return D

def getLocations(reactionData):
	D = {}
	for reaction in reactionData:
		complexId = reaction["id"]
		for element in reaction["elements"]:
			if element["molecule"] == complexId:
				loc = element["location"]
				if loc == "x" or loc == "z":
					loc = "c"
				D[complexId] = loc
	return D

def getMonomerLocationsFromOurData():
	data = read_tsv(MONOMER_MASS_FILE)
	D = dict([(x["id"], x["location"][0]) for x in data])
	return D

def removeBlaclistedReactions(reactionData):
	reactionDataFiltered = []
	for reaction in reactionData:
		if reaction["id"] in REACTION_ID_BLACKLIST:
			continue
		reactionDataFiltered.append(reaction)
	return reactionDataFiltered

def getMasses(idMass, reactionData):
	def skipReaction(reaction):
		if len(reaction["elements"]) == 1:
			return True, "Reaction not balanced"
		for element in reaction["elements"]:
			if element["type"] == "proteincomplex" and element["coeff"] == 1:
				continue
			if element["molecule"] not in idMass:
				return True, "Reaction contains molecule that doesn't have its mass specified"
			if element["coeff"] > 0 and element["molecule"] != reaction["id"]:
				return True, "Reaction has more than one product"
		return False, ""

	nNotFoundPrev = np.inf
	nNotFoundDiff = np.inf
	while nNotFoundDiff > 0:
		nNotFound = 0
		missingIds = []
		for reaction in reactionData:
			complexId = reaction["id"]
			skip, reason = skipReaction(reaction)
			if skip:
				nNotFound += 1
				missingIds.append(complexId + ": " + reason)
				continue

			mass = np.zeros(11)
			for element in reaction["elements"]:
				if element["molecule"] == complexId:
					continue
				mass += ((-1. * element["coeff"]) * idMass[element["molecule"]])
			idMass[complexId] = mass
		nNotFoundDiff = nNotFoundPrev - nNotFound
		nNotFoundPrev = nNotFound

def getComplexationAndEquilibriumReactions(idMass, reactionData):
	def skipReaction(reaction):
		if len(reaction["elements"]) == 1:
			return True, "Reaction not balanced"
		for element in reaction["elements"]:
			if element["molecule"] not in idMass:
				return True, "Reaction contains molecule that doesn't have its mass specified"
		return False, ""
	def isEquilibriumReaction(reaction):
		for element in reaction["elements"]:
			if element["type"] == "metabolite":
				return True
		return False
	def fixCompartments(reaction, ourLocations):
		
		L = []
		for element in reaction["elements"]:
			D = dict(element)
			loc = D["location"]
			if loc == "x" or loc == "z":
				D["location"] = "c"
			if element["type"] == "proteinmonomer":
				D["location"] = ourLocations[element["molecule"]]
			L.append(D)
		return L

	complexationReactions = {}
	equilibriumReactions = {}
	for reaction in reactionData:
		complexId = reaction["id"]
		skip, reason = skipReaction(reaction)
		if skip:
			continue
		if isEquilibriumReaction(reaction):
			equilibriumReactions[complexId + "_RXN"] = fixCompartments(reaction, ourLocations)
		else:
			complexationReactions[complexId + "_RXN"] = fixCompartments(reaction, ourLocations)

	return complexationReactions, equilibriumReactions


idName = getNames()
idMass = {}
metaboliteMass = getMetaboliteMasses()
idMass.update(metaboliteMass)
monomerMass = getMonomerMasses()
idMass.update(monomerMass)
rnaMass = getRnaMasses()
idMass.update(rnaMass)
ourLocations = getMonomerLocationsFromOurData()

jsonData = yaml.safe_load(open(ECOCYC_DUMP, "r"))
reactionDataFiltered = removeBlaclistedReactions(jsonData["complexations"])
idLocation = getLocations(reactionDataFiltered)
getMasses(idMass, reactionDataFiltered)

complexationReactions, equilibriumReactions = getComplexationAndEquilibriumReactions(
	idMass,
	reactionDataFiltered
	)
equilibriumBindingRates = getEquilibriumBindingRates()

idList = sorted([x for x in idLocation if x in idMass])

h = open(NEW_PROTEIN_COMPLEX_FILE, "w")
h.write('"name"\t"comments"\t"mw"\t"location"\t"reactionId"\t"id"\n')

for complexId in idList:
	h.write('%s\t""\t%s\t["%s"]\t"%s_RXN"\t"%s"\n' % (
		json.dumps(idName[complexId]) if complexId in idName else '""',
		json.dumps(idMass[complexId].tolist()),
		idLocation[complexId],
		complexId,
		complexId,
		)
	)
h.close()

h = open(NEW_COMPLEXATION_REACTION_FILE, "w")
h.write('"process"\t"stoichiometry"\t"id"\t"dir"\n')

for complexReactionId in sorted(complexationReactions):
	h.write('"complexation"\t%s\t"%s"\t1\n' % (
		json.dumps(complexationReactions[complexReactionId]),
		complexReactionId,
		)
	)
h.close()

h = open(NEW_EQUILIBRIUM_REACTION_FILE, "w")
h.write('"process"\t"stoichiometry"\t"id"\t"dir"\t"forward rate"\t"reverse rate"\n')

for equilibriumReactionId in sorted(equilibriumReactions):
	h.write('"equilibrium"\t%s\t"%s"\t1\t%g\t%g\n' % (
		json.dumps(equilibriumReactions[equilibriumReactionId]),
		equilibriumReactionId,
		equilibriumBindingRates.get(equilibriumReactionId, DEFAULT_EQUILIBRIUM_BINDING_RATES)[0],
		equilibriumBindingRates.get(equilibriumReactionId, DEFAULT_EQUILIBRIUM_BINDING_RATES)[1]
		)
	)
h.close()
