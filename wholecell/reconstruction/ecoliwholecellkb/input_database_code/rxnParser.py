import re

def parseReaction(reactionStr):
	match = re.match("^\[(?P<comp>.*?)\][ ]{0,1}: (?P<stoich>.*)$", reactionStr)
	if match != None:
		globalComp = match.group("comp")
		stoich = match.group("stoich")
	else:
		globalComp = ""
		stoich = reactionStr

	match = re.match("^(?P<lefts>.*) (?P<dir><*((==)|(--))>*) (?P<rights>.*)$", stoich)
	if match == None:
		raise Exception, "Invalid stoichiometry: %s." % (stoich)

	if match.group("dir") == "==>" or match.group("dir") == "-->":
		reactionDir = 1
	elif match.group("dir") == "<==" or match.group("dir") == "<--":
		reactionDir = -1
	elif match.group("dir") == "<==>" or match.group("dir") == "<-->":
		reactionDir = 0

	stoich = []

	lefts = match.group("lefts").split(" + ")
	for componentStr in lefts:
		coeff, mol, form, comp = parseReactionComponent(componentStr, globalComp)
		stoich.append({ "coeff": -coeff, "location": comp, "molecule": mol, "form": form})

	rights = match.group("rights").split(" + ")
	for componentStr in rights:
		coeff, mol, form, comp = parseReactionComponent(componentStr, globalComp)
		stoich.append({ "coeff": coeff, "location": comp, "molecule": mol, "form": form})

	return stoich, reactionDir

def parseLeakReaction(reactionStr):
	match = re.match("^\[(?P<comp>.*?)\][ ]{0,1}: (?P<stoich>.*)$", reactionStr)
	if match != None:
		globalComp = match.group("comp")
		stoich = match.group("stoich")
	else:
		globalComp = ""
		stoich = reactionStr

	match = re.match("^(?P<lefts>.*) (?P<dir><*((==)|(--))>*)$", stoich)
	if match == None:
		raise Exception, "Invalid stoichiometry: %s." % (stoich)

	if match.group("dir") == "==>" or match.group("dir") == "-->":
		reactionDir = 1
	elif match.group("dir") == "<==" or match.group("dir") == "<--":
		raise Exception, "Invalid stoichiometry for a leak reaction: %s." % (stoich)
	elif match.group("dir") == "<==>" or match.group("dir") == "<-->":
		reactionDir = 0

	stoich = []

	lefts = match.group("lefts").split(" + ")
	for componentStr in lefts:
		coeff, mol, form, comp = parseReactionComponent(componentStr, globalComp)
		stoich.append({ "coeff": -coeff, "location": comp, "molecule": mol, "form": form})

	if len(stoich) > 1:
		raise Exception, "Leak reaction can only have one metabolite: %s." % (stoich)

	return stoich, reactionDir

def parseReactionComponent(componentStr, globalComp):
	if globalComp == "":
		tmp = re.match("^(?P<coeff>\(\d*\.*\d*\) )*(?P<mol>.+?)(?P<form>:.+)*\[(?P<comp>.+)\]$", componentStr)
		if tmp == None:
			raise Exception, "Invalid stoichiometry: %s." % (componentStr)
		if tmp.group("coeff") == None:
			coeff = 1.0
		else:
			coeff = float(tmp.group("coeff")[1:-2])

		mol = tmp.group("mol")

		if tmp.group("form") == None:
			form = "mature"
		else:
			form = tmp.group("form")[1:]

		comp = tmp.group("comp")
	else:
		tmp = re.match("^(?P<coeff>\(\d*\.*\d*\) )*(?P<mol>.+?)(?P<form>:.+)*$", componentStr)
		if tmp == None:
			raise Exception, "Invalid stoichiometry: %s." % (componentStr)
		if tmp.group("coeff") == None:
			coeff = 1.0
		else:
			coeff = float(tmp.group("coeff")[1:-2])

		mol = tmp.group("mol")

		if tmp.group("form") == None:
			form = "mature"
		else:
			form = tmp.group("form")[1:]

		comp = globalComp

	return coeff, mol, form, comp