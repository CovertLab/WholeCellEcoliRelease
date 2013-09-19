#!/usr/bin/env python
import os
import csv
from json import loads

from checkForUnmodified import getReactionEnzymesDict, getManualReactionAnnotations

# Check for uncomplexed monomers in the reaction annotated by Feist.  This is 
# functionally analogous to checkForUnmodified.py.

def getComplexComposition():
	with open(os.path.join(os.environ['PARWHOLECELLPY'], 'data', 'parsed', 'proteinComplexes.csv'), 'rb') as csvFile:
		return {row['Frame ID']:loads(row['Composition'])['reactant'] for row in csv.DictReader(csvFile, delimiter = '\t', quotechar = '"')}

def main():
	complexComposition = getComplexComposition()
	subunits = set(subunit for composition in complexComposition.values() for subunit in composition)

	subunitAssociation = {}
	for cmplx, composition in complexComposition.items():
		for subunit in composition:
			try:
				subunitAssociation[subunit].add(cmplx)

			except KeyError:
				subunitAssociation[subunit] = {cmplx}

	homopolymers = {}
	fullComplex = {}
	partialComplex = {}

	for frameID, enzymes in getReactionEnzymesDict().items():
		if enzymes & subunits:
			key = frameID.replace('FEIST_', '')

			if len(enzymes) == 1:
				enzyme, = enzymes
				homopolymers[key] = (enzyme, subunitAssociation[enzyme])

			else:
				possibleComplexes = [cmplx for enzyme in enzymes  if subunitAssociation.has_key(enzyme) for cmplx in subunitAssociation[enzyme]]

				fullComplexes = [cmplx for cmplx in possibleComplexes if enzymes.issuperset(complexComposition[cmplx])]

				if fullComplexes:
					fullComplex[key] = fullComplexes

				else:
					partialComplex[key] = possibleComplexes

	return homopolymers, fullComplex, partialComplex

# TODO: for those that are returned, determine which correspond to homopolymers,
# since those reactions could be amended largely automatically

# TODO: could also load complexation stoichiometry and use that to annotate 
# automatically

# Possible cases:
# 1) Single monomer is associated with a complex (homopolymer)
#	- Replace monomer with complex
# 2) All monomers define a specific complex
#	- Likely a misannotation in Feist ('or' instead of 'and')
#	- Replace with complex
# 3) Subset of monomers define a specific complex
#	- Likely a misannotation in Feist ('or' instead of 'and')
#	- Replace subset with complex
# 4) Monomer(s) define part of a complex
#	- Record and output