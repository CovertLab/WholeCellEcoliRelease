from __future__ import absolute_import, division, print_function

import os
import re
import yaml
import json
from typing import Dict, List

from reconstruction.spreadsheets import read_tsv
import six


REACTIONS_FILE = os.path.join("reconstruction", "ecoli", "flat", "reactions.tsv")
ECOCYC_DUMP = os.path.join("reconstruction", "ecoli", "flat", "ecocyc_20_1_gem.json")

NEW_REACTIONS_FILE = os.path.join("reconstruction", "ecoli", "flat", "reactions_new.tsv")



def getReactionIds():
	data = read_tsv(REACTIONS_FILE)
	L = sorted(x["reaction id"] for x in data)
	return L

def getReactionStoich():
	data = read_tsv(REACTIONS_FILE)
	D = dict((x["reaction id"], x["stoichiometry"]) for x in data)
	return D

def getReactionReversibility():
	data = read_tsv(REACTIONS_FILE)
	D = dict((x["reaction id"], x["is reversible"]) for x in data)
	return D

def truncateNameRxnTrail(name):
	if name.startswith("TRANS-RXN"):
		return None
	return name.split("-RXN")[0] + "-RXN"

def truncateNameRxnStart(name):
	if name.startswith("TRANS-RXN"):
		return None
	R = re.search("(RXN0?-[0-9]*)[\-\[]", name)
	if R is not None:
		return R.groups()[0]
	return None

def truncateNameRxnTrans(name):
	R = re.search("(TRANS-RXN0?-[0-9]*)[\-\[A-Z]", name)
	if R is not None:
		return R.groups()[0]
	return None


def addFilteredEntries(rxnNamesEnzymes):
	# type: (Dict[str, List[str]]) -> None
	d = {}  # type: Dict[str, List[str]]
	def add_enzymes(name):
		# type: (str) -> None
		if name is not None:
			d[name] = enzymes

	for rxnName, enzymes in six.viewitems(rxnNamesEnzymes):
		if rxnName.endswith("-RXN"):
			continue
		add_enzymes(truncateNameRxnTrail(rxnName))
		add_enzymes(truncateNameRxnStart(rxnName))
		add_enzymes(truncateNameRxnTrans(rxnName))

	rxnNamesEnzymes.update(d)

reactionIds = getReactionIds()
reactionStoich = getReactionStoich()
reactionReversibility = getReactionReversibility()

jsonData = yaml.safe_load(open(ECOCYC_DUMP, "r"))

rxnNamesEnzymes = {
	x["name"]: x["annotation"]["enzymes"]
	for x in jsonData["reactions"] if "enzymes" in x["annotation"]}  # type: Dict[str, List[str]]

addFilteredEntries(rxnNamesEnzymes)

rxnNames = sorted(rxnNamesEnzymes)

reactionIdsSet = set(reactionIds)
rxnNamesSet = set(rxnNames)

reactionIdsEnzymes = {}
notIdenticalList = []
notFoundList = []
emptyEnzymeList = []

for reactionId in reactionIds:

	enzymeList = None

	if reactionId in rxnNamesEnzymes:
		enzymeList = rxnNamesEnzymes[reactionId]
		reactionIdsEnzymes[reactionId] = enzymeList

	elif truncateNameRxnTrail(reactionId) in rxnNamesEnzymes:
		enzymeList = rxnNamesEnzymes[truncateNameRxnTrail(reactionId)]
		reactionIdsEnzymes[reactionId] = enzymeList
		notIdenticalList.append(reactionId)

	elif truncateNameRxnStart(reactionId) in rxnNamesEnzymes:
		enzymeList = rxnNamesEnzymes[truncateNameRxnStart(reactionId)]
		reactionIdsEnzymes[reactionId] = enzymeList
		notIdenticalList.append(reactionId)

	elif truncateNameRxnTrans(reactionId) in rxnNamesEnzymes:
		enzymeList = rxnNamesEnzymes[truncateNameRxnTrans(reactionId)]
		reactionIdsEnzymes[reactionId] = enzymeList
		notIdenticalList.append(reactionId)

	else:
		notFoundList.append(reactionId)

	if enzymeList is not None and len(enzymeList) == 0:
		emptyEnzymeList.append(reactionId)


h = open(NEW_REACTIONS_FILE, "w")
h.write('"reaction id"\t"stoichiometry"\t"is reversible"\t"catalyzed by"\n')

for reactionId in reactionIds:
	enzymeList = []
	if reactionId in reactionIdsEnzymes:
		enzymeList = reactionIdsEnzymes[reactionId]
	try:
		h.write('"%s"\t%s\t%s\t%s\n' % (
			reactionId,
			json.dumps(reactionStoich[reactionId]),
			json.dumps(reactionReversibility[reactionId]),
			json.dumps(enzymeList),
			)
		)
	except Exception:
		import ipdb; ipdb.set_trace()
		raise
h.close()
