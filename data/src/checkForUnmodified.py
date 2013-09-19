#!/usr/bin/env python
import os
from csv import DictReader

# Check for unmodified proteins/RNAs/complexs in catalyzing the reactions
# reported by Feist.  This script is probably redundant in purpose with
# checkForModFormsInReactionRules.py.

MODIFIED_FILES = ('proteinComplexes_modified', 'rna_modified', 'proteinMonomers_modified')

def getUnmodifiedForm(fileName):
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', fileName + '.csv'), 'rb') as csvFile:
		return set(row['Unmodified Form'] for row in DictReader(csvFile, delimiter = '\t', quotechar = '"'))

def getAllUnmodified():
	return reduce(set.union, (getUnmodifiedForm(fileName) for fileName in MODIFIED_FILES))

def getReactionEnzymesDict():
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'reactions.csv'), 'rb') as csvFile:
		return {row['Frame ID']:loadCleanedJson(row['Enzyme']) for row in DictReader(csvFile, delimiter = '\t', quotechar = '"')}

def loadCleanedJson(jsonString):
	from json import loads

	jsonOutput = loads(jsonString)

	if jsonOutput is None:
		return set()

	else:
		return set(entry for entries in jsonOutput for entry in entries if entry)

def getManualReactionAnnotations():
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'interm_manual', 'reaction_enzyme_association.csv'), 'rb') as csvFile:
		return set(row['Reaction ID'] for row in DictReader(csvFile, delimiter = '\t', quotechar = '"'))

def main(verbose = True):
	unmodifiedEnzymes = getAllUnmodified()
	manuallyAnnotatedReactions = getManualReactionAnnotations()

	# Compile a list of reactions that utilize unmodified enzymes
	suspectReactions = {}

	for frameID, enzymes in getReactionEnzymesDict().items():
		if enzymes & unmodifiedEnzymes:
			key = frameID.replace('FEIST_', '')
			value = (enzymes & unmodifiedEnzymes, enzymes - unmodifiedEnzymes)

			suspectReactions[key] = value

			if verbose:
				if key in manuallyAnnotatedReactions:
					print key + ' [manually annotated]' + ': ' + ', '.join(value[0]) + ' (' + ', '.join(value[1]) + ')'

				else:
					print key + ': ' + ', '.join(value[0]) + ' (' + ', '.join(value[1]) + ')'

	allReactions = set(suspectReactions.keys())
	allEnzymes = reduce(set.union, (value[0] for value in suspectReactions.values()) )

	if verbose:
		print ', '.join(allEnzymes)

	return suspectReactions, allReactions, allEnzymes