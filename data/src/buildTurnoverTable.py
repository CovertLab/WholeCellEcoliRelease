#!/usr/bin/env python

import csv
import os
import json
from SOAPpy import WSDL
import numpy
import cPickle
import time

import socket
import SOAPpy

DEFAULT_BRENDA_CACHE = 'brendaCache'

EXTENSION = '.cPickleCache'

WSDL_URL = "http://www.brenda-enzymes.org/soap2/brenda.wsdl"

def buildEnzymeDict():
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

	return enzymeDict

def buildFromCache(fileName = DEFAULT_BRENDA_CACHE):
	# Uses cached BRENDA output to assemble the turnover table
	cachedOutput = loadCachedBrenda(fileName)

	enzymeDict = buildEnzymeDict()

	# Get kinetics that can be automatically acquired
	keys = sorted(enzymeDict.keys())

	total_entries = 0

	for idx, key in enumerate(keys):
		e = enzymeDict[key]

		for i in range(e.reactionCount):
			# Set reverse rate if known to only progress forward
			if e.direction[i] == 'forward only':
				e.reverseTurnover[i] = 0

				# if e.comments[i] is None:
				# 	e.comments[i] = ""

				# e.comments[i] += 'Forward only, reverse kinetics set to zero.'

			e.forwardTurnover[i], e.comments[i] = getTurnoverFromEntries( parseBrendaToEntries(cachedOutput['getTurnoverNumber'][e.EC[i]]) )
			e.kM[i] = getKmFromEntries( parseBrendaToEntries(cachedOutput['getKmValue'][e.EC[i]]) )

			if e.forwardTurnover[i]:
				total_entries += 1

			# TODO: assign forwardTurnoverUnits; I think BRENDA uses consistent units across values
			# TODO: build and output a report of unassigned/partially assigned/otherwise problematic output

			print '{frameId}[reaction #{reaction}]: ToN = {turnoverNumber}, has KM = {hasKM}'.format(
				frameId = e.frameId,
				reaction = i,
				turnoverNumber = e.forwardTurnover[i],
				hasKM = (e.kM[i] is not None)
				)

	print '{} entries'.format(total_entries)

	writeOutput(enzymeDict)

def writeOutput(enzymeDict):
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'turnover_annotation.csv'),'wb') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='"')

		csvwriter.writerow(['Enzyme Frame ID', 'EC', 'Reaction ID', 'Reaction stoichiometry', 'Direction', 'Forward', 'Units', 'Reverse', 'Units', 'Comments', 'Km'])

		keys = enzymeDict.keys()
		keys.sort()

		for key in keys:
			e = enzymeDict[key]

			for i in range(e.reactionCount):
				csvwriter.writerow([json.dumps(e.frameId), e.EC[i], e.reacID[i], e.reacStoich[i], e.direction[i], e.forwardTurnover[i], e.forwardTurnoverUnits[i], e.reverseTurnover[i], e.reverseTurnoverUnits[i], e.comments[i], json.dumps(e.kM[i])])

def parseBrendaToEntries(brendaOutput):
	# Parse output from a single webservice call into a list of dictionaries

	return [dict(field.split('*', 1) for field in entry.split('#')[:-1]) for entry in brendaOutput.split('!') if entry]

def filterEntries(entries, filterFunctions):
	# Filters a list of entries into categories that can be used to choose the best entry.

	# filterFunctions is a list of functions that should return True (accept) or False (reject), in order of priority
	# The returned list corresponds to the functions used to filter
	# The last list in the returned list contains the rejects

	output = [[] for i in range(len(filterFunctions) + 1)]

	for entry in entries:
		for i, function in enumerate(filterFunctions):
			if function(entry):
				output[i].append(entry)
				break

		else:
			# All rejected, add to last list
			output[-1].append(entry)

	return output

# Some simple functions for use with filterEntries(...)
isEcoli = lambda entry: 'coli' in entry['organism'].lower()
isWildtype = lambda entry: ('wild-type' in entry['commentary'].lower()) or ('wild type' in entry['commentary'].lower()) or ('wildtype' in entry['commentary'].lower())
is25dC = lambda entry: '25&deg;C' in entry['commentary'].lower()
hasValidTurnoverNumber = lambda entry: float(entry['turnoverNumber']) > 0
hasValidKm = lambda entry: float(entry['kmValue']) > 0

def getTurnoverFromEntries(entries):
	# TODO: check reaction substrate against expectation

	# TODO: choose highest amongst all organisms?
	# alternatively, provide a range of outputs corresponding to different degrees of specificity and permissivity
	filteredEntries = filterEntries(entries, [
		lambda entry: isEcoli(entry) and isWildtype(entry) and is25dC(entry) and hasValidTurnoverNumber(entry),
		lambda entry: isEcoli(entry) and is25dC(entry) and hasValidTurnoverNumber(entry),
		lambda entry: isEcoli(entry) and hasValidTurnoverNumber(entry),
		hasValidTurnoverNumber
		])

	for entryList in filteredEntries[:-1]:
		if entryList: # search for the max if there is an entry in the list
			entry = entryList[ numpy.argmax([entry['turnoverNumber'] for entry in entryList]) ]

			value = entry['turnoverNumber']
			commentary = entry['organism'] + ('; ' if (entry['organism'] and entry['commentary']) else '') + entry['commentary']

			break

	else:
		# all rejected, return empty
		value = ''
		commentary = ''

	return value, commentary

def getKmFromEntries(entries):
	# Since we're not really sure what we want to do with Km values at the moment, just
	# return a list of entries that satisfy basic conditions

	return [entry for entry in entries if isEcoli(entry) and hasValidKm(entry)]

def cacheBrenda(fileName = DEFAULT_BRENDA_CACHE, methodNames = ('getTurnoverNumber', 'getKmValue')):
	# Caches BRENDA for the appropriate webservice calls.  If it fails due to connection errors, it will attempt
	# to record everything that has been cached so far, and recover that on the next run.  Since it's unclear
	# why the webservice fails so frequently, calling this function automatically is not recommended.

	# TODO: add metadata to the cache so it can be annotated with the date an entry was acquired, and if desired
	# update old entires
	# TODO: allow adding new methods when recovering the cache (currently it should just break)

	ecNumbers = getEcNumbers()
	nEcNumbers = len(ecNumbers)

	# Try to recover the cached output
	try:
		cachedOutput = loadCachedBRENDA(fileName)

		print 'Recovered cached entries.'

	except IOError:
		# Build dictionaries of raw output
		cachedOutput = {methodName:{} for methodName in methodNames}

		print 'No cache detected.'

	client = WSDL.Proxy(WSDL_URL)

	print 'Beginning caching...'
	try:
		for i, ecNumber in enumerate(ecNumbers):
			for methodName in methodNames:
				if not cachedOutput[methodName].has_key(ecNumber):
					cachedOutput[methodName][ecNumber] = getattr(client, methodName)('ecNumber*' + ecNumber)
					print 'cached EC {} ({:.2%})'.format(ecNumber, 1.*i/nEcNumbers)

	except (socket.error, SOAPpy.HTTPError) as e:
		print 'Raised exception {} ({}), attempting to save generated output...'.format(type(e).__name__, str(e))

	else:
		print 'Finished caching, saving generated output...'

	finally:
		brendaCacheFileName = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', fileName + EXTENSION)
		
		with open(brendaCacheFileName, 'w') as cacheFile:
			cPickle.dump(cachedOutput, cacheFile, protocol = cPickle.HIGHEST_PROTOCOL)

		print 'Finished saving.'

def loadCachedBrenda(fileName = DEFAULT_BRENDA_CACHE):
	brendaCacheFileName = os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_auto', fileName + EXTENSION)

	with open(brendaCacheFileName, "r") as cacheFile:
		return cPickle.load(cacheFile)

def getEcNumbers():
	# Build a list of all relevant EC numbers to use in caching BRENDA
	ecNumbers = set()

	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed','reactions.csv')) as csvfile:
		dictreader = csv.DictReader(csvfile, delimiter='\t', quotechar='"')
		dictreader.next()

		for row in dictreader:
			if row['EC'] != '' and row['Enzyme'] != 'null':
				isozymes = json.loads(row['Enzyme'])

				ecNumbers.add(row['EC'])

	return ecNumbers