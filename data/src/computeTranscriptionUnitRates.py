#!/usr/bin/env python
from __future__ import division

import os
import json
import csv

import math
import numpy as np

from collections import defaultdict, Counter

WCM_PATH = os.environ['PARWHOLECELLPY']

PATH = {
	'interm_auto'			:os.path.join(WCM_PATH, 'data', 'interm_auto'),
	'interm_manual'			:os.path.join(WCM_PATH, 'data', 'interm_manual'),
	'raw'					:os.path.join(WCM_PATH, 'data', 'raw'),
	'parsed'				:os.path.join(WCM_PATH, 'data', 'parsed'),
	'Bernstein 2002'		:os.path.join(WCM_PATH, 'data', 'raw', 'Bernstein 2002.csv'),
	'Bernstein 2002 table 6':os.path.join(WCM_PATH, 'data', 'raw', 'Bernstein 2002 table 6.csv'),
	'Ecocyc_genes'			:os.path.join(WCM_PATH, 'data', 'raw', 'Ecocyc_genes.csv'),
	'gene_frameId_synonyms'	:os.path.join(WCM_PATH, 'data', 'interm_auto', 'gene_frameId_synonyms.json'),
	'transcriptionUnits'	:os.path.join(WCM_PATH, 'data', 'parsed', 'transcriptionUnits.csv'),
	'genes'					:os.path.join(WCM_PATH, 'data', 'parsed', 'genes.csv'),
	'Blattner 2005'			:os.path.join(WCM_PATH, 'data', 'raw', 'Blattner 2005.csv'),
	}

FRAMEID_SYNONYMS = json.load(open(PATH['gene_frameId_synonyms'], 'rb'))

BERNSTEIN_ID_COLUMN = 0
BERNSTEIN_HALFLIFE_COLUMN = 4
BERNSTEIN_START_ROW = 9

BERNSTEIN_6_ID_COLUMN = 0
BERNSTEIN_6_ABUNDANCE_COLUMN = 4
BERNSTEIN_6_START_ROW = 10

LOG2 = np.log(2)

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
	halfLives = {}

	with open(PATH['Bernstein 2002'], 'rb') as csvFile:
		reader = csv.reader(csvFile, dialect = 'excel-tab')#, delimiter = '\t', quotechar = '"')

		for i, row in enumerate(reader):
			sourceId = row[BERNSTEIN_ID_COLUMN].lower()
			value = row[BERNSTEIN_HALFLIFE_COLUMN]
			if i >= BERNSTEIN_START_ROW and FRAMEID_SYNONYMS.has_key(sourceId) and value:
				gene = FRAMEID_SYNONYMS[sourceId]
				halfLives[gene] = float(value)

	return halfLives

def getRelativeAbundances():
	abundances = {}

	with open(PATH['Bernstein 2002 table 6'], 'rb') as csvFile:
		reader = csv.reader(csvFile, dialect = 'excel-tab')

		for i, row in enumerate(reader):
			sourceId = row[BERNSTEIN_6_ID_COLUMN].lower()
			value = row[BERNSTEIN_6_ABUNDANCE_COLUMN]

			if i >= BERNSTEIN_6_START_ROW and FRAMEID_SYNONYMS.has_key(sourceId) and value:
				gene = FRAMEID_SYNONYMS[sourceId]
				abundances[gene] = 2**float(value)

	return abundances

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

def getUnitGroups():
	# Transcription unit groups ('unitGroups') are groups of transcription
	# units that encode the exact same set of genes

	transcriptionUnits = {}
	with open(PATH['transcriptionUnits'], 'rb') as csvFile:
		reader = csv.DictReader(csvFile, dialect = 'excel-tab')

		for row in reader:
			fId = row['Frame ID']
			genes = frozenset(json.loads(row['Genes']))

			transcriptionUnits[fId] = genes

	unitGroupsByGenes = defaultdict(list)
	for fId, genes in transcriptionUnits.items():
		unitGroupsByGenes[genes].append(fId)

	unitGroups = {
		frozenset(fIds):genes for genes, fIds in unitGroupsByGenes.items()
		}
	
	geneToGroups = defaultdict(set)
	for fIds, genes in unitGroups.items():
		for gene in genes:
			geneToGroups[gene].add(fIds)

	return unitGroups, geneToGroups

def getGeneCoords():
	geneCoords = {}
	geneDirections = {}
	with open(PATH['genes'], 'rb') as csvFile:
		reader = csv.DictReader(csvFile, dialect = 'excel-tab')

		for row in reader:
			fId = row['Frame ID']
			coord = int(row['Coordinate'])
			direction = row['Direction']

			geneCoords[fId] = coord
			geneDirections[fId] = direction

	return geneCoords, geneDirections

def main():
	unitGroupGenes, geneUnitGroups = getUnitGroups()
	uniqueGenes = set(
		gene for gene, uGs in geneUnitGroups.items() if len(uGs) == 1
		)

	geneAbundances = getRelativeAbundances()
	geneHalfLives = getHalfLives()

	# Assign abundances
	unitGroupAbundances = {}
	for fIds, genes in unitGroupGenes.items():
		unique = genes & uniqueGenes & geneAbundances.viewkeys()

		if unique:
			unitGroupAbundances[fIds] = np.mean(
				[geneAbundances[gene] for gene in unique]
				)

	# Assign half-lives
	unitGroupHalfLives = {}
	for fIds, genes in unitGroupGenes.items():
		unique = genes & uniqueGenes & geneHalfLives.viewkeys()

		if unique:
			unitGroupHalfLives[fIds] = np.mean(
				[geneHalfLives[gene] for gene in unique]
				)

	# Iteratively attempt to assign abundances
	didAssign = True
	while didAssign:
		didAssign = False

		for gene in geneAbundances.viewkeys():
			unitGroups = geneUnitGroups[gene]

			knowns = unitGroups & unitGroupAbundances.viewkeys()
			unknowns = unitGroups - unitGroupAbundances.viewkeys()

			if len(unknowns) == 1 and len(knowns) == len(unitGroups) - 1:
				unknown, = unknowns

				estimate = geneAbundances[gene] - sum(
					unitGroupAbundances[known] for known in knowns
					)

				if estimate > 0:
					unitGroupAbundances[unknown] = estimate
					didAssign = True

				else:
					# Change transcription unit abundances to account for discrepancy
					
					# Disitribute the error "evenly"
					error = estimate/sum(unitGroupAbundances[known] for known in knowns)

					for known in knowns:
						unitGroupAbundances[known] += error * unitGroupAbundances[known]

						assert unitGroupAbundances[known] >= 0, 'Attempted to assign a negative abundance'

					unitGroupAbundances[unknown] = 0
					didAssign = True

	# Iteratively attempt to assign half-lives
	didAssign = True
	while didAssign:
		didAssign = False

		for gene in (geneHalfLives.viewkeys() & geneAbundances.viewkeys()):
			unitGroups = geneUnitGroups[gene]

			knowns = unitGroups & unitGroupAbundances.viewkeys() & unitGroupHalfLives.viewkeys()
			unknowns = unitGroups - unitGroupHalfLives.viewkeys()

			if len(unknowns) == 1 and len(knowns) == len(unitGroups) - 1:
				unknown, = unknowns

				if not (unitGroupAbundances[unknown] > 0):
					continue

				ratios = [unitGroupAbundances[known]/unitGroupAbundances[unknown] for known in knowns]

				estimate = -geneHalfLives[gene] * np.log(
					0.5*(1+sum(ratios)) - sum(
						ratio * 0.5 ** (geneHalfLives[gene]/unitGroupHalfLives[known])
						for ratio, known in zip(ratios, knowns)
						)
					) / LOG2

				if estimate > 0:
					unitGroupHalfLives[unknown] = estimate
					didAssign = True

				else:
					# Change transcription unit abundances to account for discrepancy
					adjusted_abundance = sum(
						ratio*(1-0.5**(geneHalfLives[gene]/unitGroupHalfLives[known]))
						for ratio, known in zip(ratios, knowns)
						)

					adjusted_ratios = [unitGroupAbundances[known]/adjusted_abundance for known in knowns]

					adjusted_estimate = -geneHalfLives[gene] * np.log(
						0.5*(1+sum(adjusted_ratios)) - sum(
							ratio * 0.5 ** (geneHalfLives[gene]/unitGroupHalfLives[known])
							for ratio, known in zip(adjusted_ratios, knowns)
							)
						) / LOG2

					if adjusted_abundance <= 0 or adjusted_estimate <= 0 or np.isnan(adjusted_abundance) or np.isnan(adjusted_estimate):
						# Encountered a value that can't readily be calculated.
						# Since there's no easy way to address this, just set
						# the abundance to zero.

						unitGroupAbundances[unknown] = 0

					else:
						unitGroupAbundances[unknown] = adjusted_abundance					
						unitGroupHalfLives[unknown] = adjusted_estimate

						didAssign = True

	# For those remaining, assign the average half-life of known genes
	for unitGroup in (unitGroupGenes.viewkeys() - set(unitGroupHalfLives.keys())):
		genes = unitGroupGenes[unitGroup] & geneHalfLives.viewkeys()

		if genes:
			unitGroupHalfLives[unitGroup] = np.mean([geneHalfLives[gene] for gene in genes])

	# Otherwise, assign the average half-life of all transcription units
	average = np.mean(unitGroupHalfLives.values())

	for unitGroup in (unitGroupGenes.viewkeys() - set(unitGroupHalfLives.keys())):
		unitGroupHalfLives[unitGroup] = average

	return locals()


def main_old():
	transcriptionUnits, geneTUs, promoterTUs = getTranscriptionUnits()
	uniqueGenes = set(gene for gene, tUs in geneTUs.items() if len(tUs) == 1)
	
	halfLives = getHalfLives()
	assignedGeneHalfLives = set(halfLives.keys())

	relativeAbundances = getRelativeAbundances()
	# Abundances are reported according to some sort of DNA-to-RNA ratio, which
	# indicates to me that each predicted transcript abundance should be 
	# multiplied by the gene copy number.
	abundances = {gene:len(geneTUs[gene]) * relAbund for gene, relAbund in relativeAbundances.items()}
	#abundances = relativeAbundances
	assignedGeneAbundances = set(abundances.keys())

	assignedGenes = assignedGeneHalfLives & assignedGeneAbundances

	# For transcription units in which the half-life is known for at least one 
	# unique gene, assign the half-life of the transcription unit as the average
	assignedTUDegRates = set()

	for fId, tU in transcriptionUnits.items():
		genes = tU['genes'] & uniqueGenes & assignedGeneHalfLives

		if genes:
			tU['degRate'] = LOG2 / (sum(halfLives[gene] for gene in genes)/len(genes))
			assignedTUDegRates.add(fId)

	# Same for the abundances
	assignedTUAbundances = set()

	for fId, tU in transcriptionUnits.items():
		genes = tU['genes'] & uniqueGenes & assignedGeneAbundances

		if genes:
			tU['abundance'] = sum(abundances[gene] for gene in genes)/len(genes)
			#tU['abundance'] = min(abundances[gene] for gene in genes)
			assignedTUAbundances.add(fId)

	# Iteratively attempt to assign abundances
	didAssign = True
	while didAssign:
		didAssign = False

		for gene in assignedGeneAbundances:
			knownTUAbundances = geneTUs[gene] & assignedTUAbundances
			unknownTUAbundances = geneTUs[gene] - assignedTUAbundances

			if len(unknownTUAbundances) == 1 and len(knownTUAbundances) == len(geneTUs[gene]) - 1:
				unknownTU, = unknownTUAbundances

				estimate = abundances[gene] - sum(transcriptionUnits[knownTU]['abundance'] for knownTU in knownTUAbundances)

				if estimate < 0:
					#print 'negative abundance estimate'
					#print abundances[gene], estimate
					continue

				transcriptionUnits[unknownTU]['abundance'] = estimate
				assignedTUAbundances.add(unknownTU)

				didAssign = True

	# Iteratively attempt to assign degradation rates
	nOriginal = len(assignedTUDegRates)
	nAssigned = 0

	didAssign = True
	while didAssign:
		didAssign = False

		for gene in assignedGenes:
			knownTUDegRates = geneTUs[gene] & assignedTUAbundances & assignedTUDegRates
			unknownTUDegRates = geneTUs[gene] & assignedTUAbundances - assignedTUDegRates

			if len(geneTUs[gene]) == 2 and len(knownTUDegRates) == 1 and len(unknownTUDegRates) == 1:
				knownTU, = knownTUDegRates
				unknownTU, = unknownTUDegRates

				if REQUIRE_SAME_PROMOTER and transcriptionUnits[knownTU]['promoter'] != transcriptionUnits[unknownTU]['promoter']:
					continue

				knownHL = LOG2/transcriptionUnits[knownTU]['degRate']

				ratio = abundances[gene]/transcriptionUnits[knownTU]['abundance'] - 1

				#estimate = -1/halfLives[gene] * math.log(1-2**(-halfLives[gene]/knownHL))
				try:
					estimate = -1/halfLives[gene] * math.log(1/2 + 1/(2*ratio) - 1/ratio * 2**(-halfLives[gene]/knownHL))

				except:
					continue

				if estimate < 0:
					print estimate
					continue

				#print estimate
				transcriptionUnits[unknownTU]['degRate'] = estimate
				assignedTUDegRates.add(unknownTU)
				didAssign = True

				nAssigned += 1

		# for gene in assignedGeneHalfLives:
		# 	knownTUDegRates = geneTUs[gene] & assignedTUDegRates
		# 	unknownTUDegRates = geneTUs[gene] - assignedTUDegRates

		# 	if len(geneTUs[gene]) == 2 and len(knownTUDegRates) == 1 and len(unknownTUDegRates) == 1:
		# 		knownTU, = knownTUDegRates
		# 		unknownTU, = unknownTUDegRates

		# 		if REQUIRE_SAME_PROMOTER and transcriptionUnits[knownTU]['promoter'] != transcriptionUnits[unknownTU]['promoter']:
		# 			continue

		# 		knownHL = LOG2/transcriptionUnits[knownTU]['degRate']

		# 		ratio = 1

		# 		#estimate = -1/halfLives[gene] * math.log(1-2**(-halfLives[gene]/knownHL))
		# 		try:
		# 			estimate = -1/halfLives[gene] * math.log(1/2 + 1/(2*ratio) - 1/ratio * 2**(-halfLives[gene]/knownHL))

		# 		except:
		# 			continue

		# 		if estimate < 0:
		# 			print estimate
		# 			continue

		# 		#print estimate
		# 		transcriptionUnits[unknownTU]['degRate'] = estimate
		# 		assignedTUDegRates.add(unknownTU)
		# 		didAssign = True

		# 		nAssigned += 1

	print '{} originally assigned + {} new = {} total assignments'.format(nOriginal, nAssigned, nOriginal + nAssigned)

	# TODO: write rates to file
	# TODO: assign rates to unassignable

def redundant():
	# find transcription units with the same genes
	transcriptionUnits, _, _ = getTranscriptionUnits()

	genes = [frozenset(tU['genes']) for tU in transcriptionUnits.values()]

	counter = Counter(genes)

	import ipdb
	ipdb.set_trace()