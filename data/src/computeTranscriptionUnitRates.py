#!/usr/bin/env python
from __future__ import division

import os
import json
import csv

import math

from collections import defaultdict

WCM_PATH = os.environ['PARWHOLECELLPY']

PATH = {
	'interm_auto'			:os.path.join(WCM_PATH, 'data', 'interm_auto'),
	'interm_manual'			:os.path.join(WCM_PATH, 'data', 'interm_manual'),
	'raw'					:os.path.join(WCM_PATH, 'data', 'raw'),
	'parsed'				:os.path.join(WCM_PATH, 'data', 'parsed'),
	'Bernstein 2002'		:os.path.join(WCM_PATH, 'data', 'raw', 'Bernstein 2002.csv'),
	'Ecocyc_genes'			:os.path.join(WCM_PATH, 'data', 'raw', 'Ecocyc_genes.csv'),
	'gene_frameId_synonyms'	:os.path.join(WCM_PATH, 'data', 'interm_auto', 'gene_frameId_synonyms.json'),
	'transcriptionUnits'	:os.path.join(WCM_PATH, 'data', 'parsed', 'transcriptionUnits.csv'),
	'Blattner 2005'			:os.path.join(WCM_PATH, 'data', 'raw', 'Blattner 2005.csv'),
	}

FRAMEID_SYNONYMS = json.load(open(PATH['gene_frameId_synonyms'], 'rb'))

BERNSTEIN_ID_COLUMN = 0
BERNSTEIN_HALFLIFE_COLUMN = 4
BERNSTEIN_START_ROW = 9

LOG2 = math.log(2)

REQUIRE_SAME_PROMOTER = False

def getTranscriptionUnits():
	transcriptionUnits = {}

	genes = defaultdict(set)
	promoters = defaultdict(set)

	with open(PATH['transcriptionUnits'], 'rb') as csvFile:
		reader = csv.DictReader(csvFile, dialect = 'excel-tab')#, delimiter = '\t', quotechar = '"')

		for row in reader:
			fId = row['Frame ID']
			gs = set(json.loads(row['Genes']))
			p = row['Promoter']

			transcriptionUnits[fId] = {'genes':gs, 'promoter':p}

			for g in gs:
				genes[g].add(fId)

			promoters[p].add(fId)

	return transcriptionUnits, genes, promoters

def getHalfLives():
	halfLives = defaultdict(lambda: None)

	with open(PATH['Bernstein 2002'], 'rb') as csvFile:
		reader = csv.reader(csvFile, dialect = 'excel-tab')#, delimiter = '\t', quotechar = '"')

		for i, row in enumerate(reader):
			sourceId = row[BERNSTEIN_ID_COLUMN].lower()
			value = row[BERNSTEIN_HALFLIFE_COLUMN]
			if i >= BERNSTEIN_START_ROW and FRAMEID_SYNONYMS.has_key(sourceId) and value:
				gene = FRAMEID_SYNONYMS[sourceId]
				halfLives[gene] = float(value)

	return halfLives

# def getExpressionRates():
# 	expRates = defaultdict(lambda: None)

# 	with open(PATH['Blattner 2005'], 'rb') as csvFile:
# 		reader = csv.reader(csvFile, dialect = 'excel-tab')#, delimiter = '\t', quotechar = '"')

# 		for i, row in enumerate(reader):
# 			sourceId = row[1].lower()
# 			if i > 97 and FRAMEID_SYNONYMS.has_key(sourceId):
# 				gene = FRAMEID_SYNONYMS[sourceId]
# 				expRates[gene] = sum(float(value) for value in row[2:7])/5 # TODO: take the base-2 log of these values?

# 	return expRates

def main():
	transcriptionUnits, geneTUs, promoterTUs = getTranscriptionUnits()
	uniqueGenes = set(gene for gene, tUs in geneTUs.items() if len(tUs) == 1)
	
	halfLives = getHalfLives()
	# avgHalfLife = sum(halfLives.values())/len(halfLives)
	assignedGeneHalfLives = set(halfLives.keys())

	# For transcription units in which the half-life is known for at least one 
	# unique gene, assign the half-life of the transcription unit as the average
	assignedTUDegRates = set()

	for fId, tU in transcriptionUnits.items():
		genes = tU['genes'] & uniqueGenes & assignedGeneHalfLives

		if genes:
			tU['degRate'] = LOG2 / (sum(halfLives[gene] for gene in genes)/len(genes))
			assignedTUDegRates.add(fId)

	nOriginal = len(assignedTUDegRates)
	nAssigned = 0

	didAssign = True
	while didAssign:
		didAssign = False

		for gene in assignedGeneHalfLives:
			knownTUDegRates = geneTUs[gene] & assignedTUDegRates
			unknownTUDegRates = geneTUs[gene] - assignedTUDegRates

			if len(geneTUs[gene]) == 2 and len(knownTUDegRates) == 1 and len(unknownTUDegRates) == 1:
				knownTU, = knownTUDegRates
				unknownTU, = unknownTUDegRates

				if REQUIRE_SAME_PROMOTER and transcriptionUnits[knownTU]['promoter'] != transcriptionUnits[unknownTU]['promoter']:
					continue

				knownHL = LOG2/transcriptionUnits[knownTU]['degRate']

				estimate = -1/halfLives[gene] * math.log(1-2**(-halfLives[gene]/knownHL))

				#print estimate
				transcriptionUnits[unknownTU]['degRate'] = estimate
				assignedTUDegRates.add(unknownTU)
				didAssign = True

				nAssigned += 1

	print '{} originally assigned + {} new = {} total assignments'.format(nOriginal, nAssigned, nOriginal + nAssigned)

	import ipdb
	ipdb.set_trace()

	# TODO: write rates to file
	# TODO: assign rates to unassignable