#!/usr/bin/env python
from __future__ import division

import os
import json
import csv

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

def getBlattnerAbundances():
	abundances = {}

	with open(PATH['Blattner 2005'], 'rb') as csvFile:
		reader = csv.reader(csvFile, dialect = 'excel-tab')

		for i, row in enumerate(reader):
			sourceId = row[1].lower()
			if i > 97 and FRAMEID_SYNONYMS.has_key(sourceId):
				gene = FRAMEID_SYNONYMS[sourceId]
				abundances[gene] = np.sum([float(value) for value in row[2:7]])/5

	return abundances

def calculateHalfLives(unitGroupGenes, geneUnitGroups, uniqueGenes):
	geneAbundances = getRelativeAbundances()
	geneHalfLives = getHalfLives()

	# Records of transcription unit group assignment success
	recordSimple = set() # calculatable from unique genes
	recordDeconvolved = set() # successfully determined by comparing with related transcription units
	recordContradictory = set() # required serious adjustment
	recordAmbiguous = set() # couldn't separate, not enough information
	recordUnknown = set() # literally no genes with known values
	records = {
		'simple':recordSimple,
		'deconvolved':recordDeconvolved,
		'contradictory':recordContradictory,
		'ambiguous':recordAmbiguous,
		'unknown':recordUnknown,
		}

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

	recordSimple |= set(unitGroupHalfLives.viewkeys())

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

				estimate = geneAbundances[gene] - np.sum(
					[unitGroupAbundances[known] for known in knowns]
					)

				if estimate > 0:
					unitGroupAbundances[unknown] = estimate
					didAssign = True

					recordDeconvolved.add(unknown)

				else:
					# Change transcription unit abundances to account for discrepancy
					
					# Disitribute the error "evenly"
					error = estimate/np.sum([unitGroupAbundances[known] for known in knowns])

					for known in knowns:
						unitGroupAbundances[known] += error * unitGroupAbundances[known]

						assert unitGroupAbundances[known] >= 0, 'Attempted to assign a negative abundance'

					unitGroupAbundances[unknown] = 0
					didAssign = True

					recordContradictory.add(unknown)
					recordContradictory |= knowns

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
					0.5*(1+np.sum(ratios)) - np.sum(
						[ratio * 0.5 ** (geneHalfLives[gene]/unitGroupHalfLives[known])
						for ratio, known in zip(ratios, knowns)]
						)
					) / LOG2

				if estimate > 0:
					unitGroupHalfLives[unknown] = estimate
					didAssign = True

					recordDeconvolved.add(unknown)

				else:
					# Change transcription unit abundances to account for discrepancy
					adjusted_abundance = np.sum(
						[ratio*(1-0.5**(geneHalfLives[gene]/unitGroupHalfLives[known]))
						for ratio, known in zip(ratios, knowns)]
						)

					adjusted_ratios = [unitGroupAbundances[known]/adjusted_abundance for known in knowns]

					adjusted_estimate = -geneHalfLives[gene] * np.log(
						0.5*(1+np.sum(adjusted_ratios)) - np.sum(
							[ratio * 0.5 ** (geneHalfLives[gene]/unitGroupHalfLives[known])
							for ratio, known in zip(adjusted_ratios, knowns)]
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

					recordContradictory.add(unknown)
					recordContradictory |= knowns

	recordSimple -= recordContradictory
	recordDeconvolved -= recordContradictory

	# For those remaining, assign the average half-life of known genes
	for unitGroup in (unitGroupGenes.viewkeys() - set(unitGroupHalfLives.keys())):
		genes = unitGroupGenes[unitGroup] & geneHalfLives.viewkeys()

		if genes:
			unitGroupHalfLives[unitGroup] = np.mean([geneHalfLives[gene] for gene in genes])

			recordAmbiguous.add(unitGroup)

	# Otherwise, assign the average half-life of all transcription units
	average = np.mean(unitGroupHalfLives.values())

	for unitGroup in (unitGroupGenes.viewkeys() - set(unitGroupHalfLives.keys())):
		unitGroupHalfLives[unitGroup] = average

		recordUnknown.add(unitGroup)

	recordDeconvolved -= recordAmbiguous
	recordDeconvolved -= recordUnknown

	recordContradictory -= recordAmbiguous
	recordContradictory -= recordUnknown

	return unitGroupHalfLives, records

def calculateExpressionRates(unitGroupGenes, geneUnitGroups, uniqueGenes, unitGroupHalfLives):
	geneAbundances = getBlattnerAbundances()

	unitGroupAbundances = {}
	for unitGroup, genes in unitGroupGenes.items():
		unique = genes & uniqueGenes & geneAbundances.viewkeys()

		if unique:
			unitGroupAbundances[unitGroup] = np.sum(
				[geneAbundances[gene] for gene in unique]
				)/len(unique)

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

				estimate = geneAbundances[gene] - np.sum(
					[unitGroupAbundances[known] for known in knowns]
					)

				if estimate > 0:
					unitGroupAbundances[unknown] = estimate
					didAssign = True

				else:
					# Change transcription unit abundances to account for discrepancy
					
					# Disitribute the error "evenly"
					error = estimate/np.sum([unitGroupAbundances[known] for known in knowns])

					for known in knowns:
						unitGroupAbundances[known] += error * unitGroupAbundances[known]

						assert unitGroupAbundances[known] >= 0, 'Attempted to assign a negative abundance'

					unitGroupAbundances[unknown] = 0
					didAssign = True

	# For those remaining, assign the average abundance of known genes
	for unitGroup in (unitGroupGenes.viewkeys() - set(unitGroupAbundances.keys())):
		genes = unitGroupGenes[unitGroup] & geneAbundances.viewkeys()

		if genes:
			unitGroupAbundances[unitGroup] = np.mean([geneAbundances[gene] for gene in genes])

	# Otherwise, assign the average abundance of all transcription units
	average = np.mean(unitGroupAbundances.values())

	for unitGroup in (unitGroupGenes.viewkeys() - set(unitGroupAbundances.keys())):
		unitGroupAbundances[unitGroup] = average

	# Use the abundances and the half-lives to determine expression rates
	unitGroupExpressionRates = {}
	for unitGroup in unitGroupGenes:
		unitGroupExpressionRates[unitGroup] = unitGroupAbundances[unitGroup] * LOG2 / unitGroupHalfLives[unitGroup]

	return unitGroupExpressionRates

def write(unitGroupHalfLives, unitGroupExpressionRates):
	rows = []
	with open(PATH['transcriptionUnits'], 'rb') as csvFile:
		reader = csv.reader(csvFile, dialect = 'excel-tab')

		for row in reader:
			rows.append(row)

	for i, row in enumerate(rows):
		if i == 0:
			row.append('Degradation rate (1/min)')
			row.append('Expression rate (a.u./min)')

		else:
			fId = row[0]
			row.append(LOG2/unitGroupHalfLives[fId])
			row.append(unitGroupExpressionRates[fId])

	with open(os.path.join(WCM_PATH, 'data', 'parsed', 'transcriptionUnits_with_rates.csv'), 'wb') as csvFile:
		writer = csv.writer(csvFile, dialect = 'excel-tab')

		for row in rows:
			writer.writerow(row)

def unitGroupsToTranscriptionUnits(dct):
	return {
		transcriptionUnit:value
		for unitGroup, value in dct.items()
		for transcriptionUnit in unitGroup
		}

def main():
	unitGroupGenes, geneUnitGroups = getUnitGroups()
	uniqueGenes = set(
		gene for gene, uGs in geneUnitGroups.items() if len(uGs) == 1
		)

	unitGroupHalfLives, unused = calculateHalfLives(unitGroupGenes, geneUnitGroups,
		uniqueGenes)

	unitGroupExpressionRates = calculateExpressionRates(unitGroupGenes,
		geneUnitGroups, uniqueGenes, unitGroupHalfLives)

	write(
		unitGroupsToTranscriptionUnits(unitGroupHalfLives),
		unitGroupsToTranscriptionUnits(unitGroupExpressionRates)
		)

	#return Bunch(**locals())

from textwrap import fill as tw_fill

def makeReport():
	unitGroupGenes, geneUnitGroups = getUnitGroups()
	uniqueGenes = set(
		gene for gene, uGs in geneUnitGroups.items() if len(uGs) == 1
		)

	unitGroupHalfLives, records = calculateHalfLives(unitGroupGenes, geneUnitGroups,
		uniqueGenes)

	output = ''

	output += '''Note that all counts are on a per-transcription-unit-group basis, where a group
is defined as all transcription units with the same set of genes.

Summary statistics
{} simply assigned (Cat. 1)
{} were converted from gene half-lives to transcription unit half-lives (Cat. 2)
{} were contradictory (half-lives did not agree with TU structures) (Cat. 3)
{} were assigned based on the average of known gene half-lives (Cat. 4)
{} had no known gene half-lives and were assigned according to the average (Cat. 5)'''.format(
		len(records['simple']), len(records['deconvolved']),
		len(records['contradictory']), len(records['ambiguous']),
		len(records['unknown'])
		)

	for n, key in enumerate(['simple', 'deconvolved', 'contradictory', 'ambiguous', 'unknown']):
		output += '\n\n'

		output += 'Category {}:\n'.format(n+1)
		output += '#'*79 + '\n'

		output += '\n'.join(
			tw_fill(
				'{}: ({})'.format(', '.join(tU for tU in tUs), ', '.join(gene for gene in unitGroupGenes[tUs])),
				width = 79,
				subsequent_indent = '\t',
				)
			for tUs in records[key]
			)

	with open(os.path.join(WCM_PATH, 'data', 'parsed', 'report_for_Nick.txt'), 'wb') as report:
		report.write(output)

	#print output


class Bunch(object):
	'''
	Helper class for development/debugging.  Usage:
	return Bunch(**locals())
	'''
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)

if __name__ == '__main__':
	main()