#!/usr/bin/env python

import csv
import os
import json
from SOAPpy import WSDL


def createTurnoverTable():
	enzymeDict = {}
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed','reactions.csv')) as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		dictreader.next()
		for row in dictreader:
			if row['EC'] != '' and row['Enzyme'] != 'null':
				isozymes = json.loads(row['Enzyme'])
				for iso in isozymes:
					iso = tuple(iso)
					if not enzymeDict.has_key(iso):
						newEnzyme = enzyme()
						enzymeDict[iso] = newEnzyme
					enzymeDict[iso].frameId = iso
					enzymeDict[iso].EC.append(row['EC'])
					enzymeDict[iso].reacID.append(row['Frame ID'])
					enzymeDict[iso].reacStoich.append(row['Stoichiometry (pH 7.2)'])
					enzymeDict[iso].direction.append(row['Direction'])

	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'intermediate', 'turnover_annotation.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['Enzyme Frame ID', 'EC', 'Reaction ID', 'Reaction stoichiometry', 'Direction', 'Turnover Forward (s^-1)', 'Turnover Reverse (s^-1)', 'Comments'])

		keys = enzymeDict.keys()
		keys.sort()

		for key in keys:
			e = enzymeDict[key]
			for i in range(len(e.EC)):
				reverseSet = None
				if e.direction[i] == 'forward only':
					reverseSet = 0
				csvwriter.writerow([json.dumps(e.frameId), e.EC[i], e.reacID[i], e.reacStoich[i], e.direction[i], None, reverseSet])

class enzyme():
	def __init__(self):
		self.frameId = None
		self.EC = []
		self.reacID = []
		self.reacStoich = []
		self.direction = []








def getBRENDA():
	# Build complete list of EC numbers in Feist
	ECnumbers = []
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'raw','Feist_reactions.csv')) as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')

		for row in dictreader:
			if row['proteinClass'] != '':
				ECnumbers.append(row['proteinClass'])
		ECnumbers = set(ECnumbers)
		ECnumbers = list(ECnumbers)

	# Query BRENDA for data on all EC numbers
	reac = parseBrendaReaction("ecNumber*2.3.1.40#organism*Escherichia coli")
	result = client.getTurnoverNumber("ecNumber*2.3.1.40#organism*Escherichia coli")
	ipdb.set_trace()

def parseBrendaReaction(line):
	reac = BRENDA_reaction()
	return parseBrendaEntry(reac,line)

def parseBrendaTurnover(line):
	turn = BRENDA_turnover()
	return parseBrenda(turn,line)

def parseBrendaEntry(obj, line):
	wsdl = "http://www.brenda-enzymes.org/soap2/brenda.wsdl"
	client = WSDL.Proxy(wsdl)
	result = client.getReaction(line)
	entries = result.split('!')
	for e in entries:
		fields = e.split('#')
		for f in fields:
			if len(f.split('*')) > 1:
				attr_name = f.split('*')[0]
				attr_value = f.split('*')[1]
				setattr(obj, attr_name, attr_value)
	return obj

class BRENDA_reaction():
	def __init__(self):
		self.ecNumber = None
		self.reaction = None
		self.commentary = None
		self.literature = None
		self.organism = None

class BRENDA_turnover():
	def __init__(self):
		pass
