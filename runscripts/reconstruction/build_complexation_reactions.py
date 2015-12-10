import os
from reconstruction.spreadsheets import JsonReader
import csv
import yaml
import json
import numpy as np

CSV_DIALECT = csv.excel_tab

METABOLITE_MASS_FILE = os.path.join("reconstruction", "ecoli", "flat", "metabolites.tsv")
MONOMER_MASS_FILE = os.path.join("reconstruction", "ecoli", "flat", "proteins.tsv")
RNA_MASS_FILE = os.path.join("reconstruction", "ecoli", "flat", "rnas.tsv")
WATER_MASS_FILE = os.path.join("reconstruction", "ecoli", "flat", "water.tsv")
EXISTING_PROTEIN_COMPLEX_FILE = os.path.join("reconstruction", "ecoli", "flat", "proteinComplexes.tsv")
NEW_PROTEIN_COMPLEX_FILE = os.path.join("reconstruction", "ecoli", "flat", "proteinComplexes_new.tsv")

ECOCYC_DUMP = os.path.join("reconstruction", "ecoli", "flat", "eco_wc_test_fun_truncated.json")




def getMetaboliteMasses():
	reader = JsonReader(open(METABOLITE_MASS_FILE, "r"), dialect = CSV_DIALECT)
	data = [row for row in reader]
	D = dict([(x["id"].encode("utf-8"), np.array([0, 0, 0, 0, 0, 0, 0, x["mw7.2"], 0, 0, 0])) for x in data])
	reader = JsonReader(open(WATER_MASS_FILE, "r"), dialect = CSV_DIALECT)
	data = [row for row in reader]
	E = dict([(x["id"].encode("utf-8"), np.array([0, 0, 0, 0, 0, 0, 0, 0, x["mw7.2"], 0, 0])) for x in data])
	return dict(D, **E)

def getMonomerMasses():
	reader = JsonReader(open(MONOMER_MASS_FILE, "r"), dialect = CSV_DIALECT)
	data = [row for row in reader]
	D = dict([(x["id"].encode("utf-8"), np.array(x["mw"])) for x in data])
	return D

def getRnaMasses():
	reader = JsonReader(open(RNA_MASS_FILE, "r"), dialect = CSV_DIALECT)
	data = [row for row in reader]
	D = dict([(x["id"].encode("utf-8"), np.array(x["mw"])) for x in data])
	return D

def getNames():
	reader = JsonReader(open(EXISTING_PROTEIN_COMPLEX_FILE, "r"), dialect = CSV_DIALECT)
	data = [row for row in reader]
	D = dict([(x["id"].encode("utf-8"), x["name"].encode("utf-8")) for x in data])
	return D

def getLocations(reactionData):
	D = {}
	for reaction in reactionData:
		complexId = reaction["id"]
		for element in reaction["elements"]:
			if element["molecule"] == complexId:
				D[complexId] = element["location"]
	return D

def getMasses(idMass, reactionData):
	def skipReaction(reaction):
		if len(reaction["elements"]) == 1:
			return True, "Reaction not balanced"
		for element in reaction["elements"]:
			if element["type"] == "proteincomplex" and element["coeff"] == 1:
				continue
			# if element["type"] != "proteinmonomer":
			# 	return True, "Reaction contains non-monomer element"
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
				# print "Skipping complex %s: %s" % (complexId, reason)
				nNotFound += 1
				missingIds.append(complexId + ": " + reason)
				continue

			mass = np.zeros(11)
			for element in reaction["elements"]:
				if element["molecule"] == complexId:
					continue
				mass += ((-1. * element["coeff"]) * idMass[element["molecule"]])
			idMass[complexId] = mass
		print nNotFound
		nNotFoundDiff = nNotFoundPrev - nNotFound
		nNotFoundPrev = nNotFound

idName = getNames()
idMass = {}
metaboliteMass = getMetaboliteMasses()
idMass.update(metaboliteMass)
monomerMass = getMonomerMasses()
idMass.update(monomerMass)
rnaMass = getRnaMasses()
idMass.update(rnaMass)

jsonData = yaml.load(open(ECOCYC_DUMP, "r"))
idLocation = getLocations(jsonData["complexations"])
getMasses(idMass, jsonData["complexations"])

idList = sorted([x for x in idLocation if x in idMass])

print [x for x in idList if x not in idName]

h = open(NEW_PROTEIN_COMPLEX_FILE, "w")
h.write('"name"\t"comments"\t"mw"\t"location"\t"reactionId"\t"id"\n')

for complexId in idList:
	h.write('"%s"\t""\t%s\t["%s"]\t"%s_RXN"\t"%s"\n' % (
		idName[complexId] if complexId in idName else "",
		json.dumps(idMass[complexId].tolist()),
		idLocation[complexId],
		complexId,
		complexId,
		)
	)
h.close()

