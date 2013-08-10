#!/usr/bin/env python

import csv
import os
import json
from SOAPpy import WSDL
import numpy
import cPickle
import time

def persistantRun():
	# If brenda stops responding, just keep re-running buildTurnoverTable()
	# We cache intermediate steps until it completes
	val = None
	while val != 1:
		try:
			val = buildTurnoverTable()
		except:
			time.sleep(60)
	print 'COMPLETED!'

def buildTurnoverTable(clearCache = False):
	enzymeDict = {}
	# Build all enzymes
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
					enzymeDict[iso].reactionCount += 1
					enzymeDict[iso].forwardTurnover.append(None)
					enzymeDict[iso].forwardTurnoverUnits.append(None)
					enzymeDict[iso].reverseTurnover.append(None)
					enzymeDict[iso].reverseTurnoverUnits.append(None)
					enzymeDict[iso].kM.append(None)
					enzymeDict[iso].comments.append('')

	wsdl = "http://www.brenda-enzymes.org/soap2/brenda.wsdl"
	client = WSDL.Proxy(wsdl)
	cacheFileName = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto','turnoverCache.cPickle')
	if not os.path.exists(cacheFileName) or clearCache:
		startIdx = 0
	else:
		enzymeDict, startIdx = cPickle.load(open(cacheFileName, "r"))
	cacheCount = 0
	# Get kinetics that can be automatically gotten
	keys = sorted(enzymeDict.keys())
	for idx in xrange(startIdx, len(keys)):
		e = enzymeDict[keys[idx]]
		for i in range(e.reactionCount):
			# Set reverse rate if known to only progress forward
			if e.direction[i] == 'forward only':
				e.reverseTurnover[i] = 0
				if e.comments[i] == None:
					e.comments[i] = ""
				e.comments[i] += 'Forward only, reverse kinetics set to zero.'

			e.forwardTurnover[i], e.comments[i] = parseBrendaTurnover(client, "ecNumber*" + e.EC[i])
			e.kM[i] = parseBrendaKm(client, "ecNumber*" + e.EC[i])
			print "%s[%d]: %s" % (e.frameId, i, str(e.forwardTurnover[i]))
			if e.kM[i] != None:
				print "%s[%d]: %s" % (e.frameId, i, str(e.kM[i][0]))
			cacheCount += 1
			if cacheCount % 10 == 0:
				cPickle.dump((enzymeDict, idx), open(cacheFileName, "w"), protocol = cPickle.HIGHEST_PROTOCOL)

	cPickle.dump((enzymeDict, idx), open(cacheFileName, "w"), protocol = cPickle.HIGHEST_PROTOCOL)
	# Write output
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'turnover_annotation.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['Enzyme Frame ID', 'EC', 'Reaction ID', 'Reaction stoichiometry', 'Direction', 'Forward', 'Units', 'Reverse', 'Units', 'Comments', 'Km'])

		keys = enzymeDict.keys()
		keys.sort()

		for key in keys:
			e = enzymeDict[key]
			for i in range(e.reactionCount):
				csvwriter.writerow([json.dumps(e.frameId), e.EC[i], e.reacID[i], e.reacStoich[i], e.direction[i], e.forwardTurnover[i], e.forwardTurnoverUnits[i], e.reverseTurnover[i], e.reverseTurnoverUnits[i], e.comments[i], json.dumps(e.kM[i])])
	return 1

class enzyme():
	def __init__(self):
		self.frameId = None
		self.EC = []
		self.reacID = []
		self.reacStoich = []
		self.direction = []
		self.forwardTurnover = []
		self.reverseTurnover = []
		self.forwardTurnoverUnits = []
		self.reverseTurnoverUnits = []
		self.kM = []
		self.reactionCount = 0
		self.comments = []


# def parseBrendaReaction(line):
# 	reac = BRENDA_reaction()
# 	return parseBrendaEntry(reac,line)

def parseBrendaTurnover(client, line):


	result = client.getTurnoverNumber(line)
	if len(result) == 0:
		return None, None
	entries = result.split('!')
	L = []
	for e in entries:
		L.append(dict([x.split("*", 1) for x in e.split("#") if len(x) > 0]))
	if any("coli" in x["organism"].lower() for x in L):

		for entry in L:
			possValue = []
			if entry['commentary'].lower().count('wild') and entry['commentary'].count('25'):
				possValue.append(entry['turnoverNumber'])
				ipdb.set_trace()
		
		if len(possValue):
			value = max(possValue)
		else:
			value = -1


		return (value, "In E. coli")
	maxVal, maxIdx = numpy.max([float(x["turnoverNumber"]) for x in L]), numpy.argmax([float(x["turnoverNumber"]) for x in L])
	return (maxVal, "organism: %s\tcomments:%s" % (L[maxIdx]["organism"], L[maxIdx]["commentary"]))


def parseBrendaKm(client, line):
	result = client.getKmValue(line)

	if len(result) == 0:
		return None
	entries = result.split('!')
	L = []
	for e in entries:
		L.append(dict([x.split("*", 1) for x in e.split("#") if len(x) > 0]))
	if any("coli" in x["organism"].lower() for x in L):
		return [dict([("kmValue", x["kmValue"]), ("kmValueMaximum", x["kmValueMaximum"]), ("substrate", x["substrate"]), ("organism", x["organism"]), ("commentary", x["commentary"]), ("brendaLiteratureId", x["literature"])]) for x in L if "coli" in x["organism"].lower()]
	return None


# class BRENDA_reaction():
# 	def __init__(self):
# 		self.ecNumber = None
# 		self.reaction = None
# 		self.commentary = None
# 		self.literature = None
# 		self.organism = None

class BRENDA_turnover():
	def __init__(self):
		pass
