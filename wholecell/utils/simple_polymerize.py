#!/usr/bin/env python

"""
simple_polymerize.py

A generalization of the algorithm used in our simplified (bulk) transcription
and translation submodels.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/22/2014
"""

import numpy as np

POLYMER_CREATION_FALLOFF_RATE = 0.9

def simplePolymerize(templateMonomerCounts, enzymaticLimitation,
		monomerCounts, synthesisProbabilities, randStream):

	assert np.abs(synthesisProbabilities.sum() - 1) < 1e-5

	templateLengths = templateMonomerCounts.sum(1)

	avgTemplateLength = templateLengths.mean()

	nPolymersToCreate = int(
		min(enzymaticLimitation, monomerCounts.sum()) / avgTemplateLength
		)
	
	polymersCreated = np.zeros_like(templateLengths)

	while True:
		if nPolymersToCreate <= 0:
			break

		# Can only polymerize if 1) there are enough monomers and 2) there is enough enzymatic capacity
		canPolymerize = ((monomerCounts >= templateMonomerCounts).all(1) &
			(enzymaticLimitation >= templateLengths))

		if not np.any(canPolymerize):
			break

		adjustedSynthProb = canPolymerize * synthesisProbabilities
		adjustedSynthProb /= adjustedSynthProb.sum()

		# Choose numbers of each polymer to synthesize
		nNewPolymers = randStream.mnrnd(nPolymersToCreate, adjustedSynthProb)

		# TODO: determine why this assertion fails on Hudson
		# assert nNewPolymers.sum() == nPolymersToCreate

		monomersUsed = np.dot(nNewPolymers, templateMonomerCounts)

		enzymaticCapacityUsed = np.dot(nNewPolymers, templateLengths).sum()

		# Reduce number of polymers to create if not enough monomers
		if (monomerCounts < monomersUsed).any():
			nPolymersToCreate = max(int(nPolymersToCreate * POLYMER_CREATION_FALLOFF_RATE), 1)
			continue

		# Reduce number of polymers to create if not enough enzymatic capacity
		if enzymaticLimitation < enzymaticCapacityUsed:
			nPolymersToCreate = max(int(nPolymersToCreate * POLYMER_CREATION_FALLOFF_RATE), 1)
			continue

		# Use resources and create the polymers

		enzymaticLimitation -= enzymaticCapacityUsed

		monomerCounts -= monomersUsed

		polymersCreated += nNewPolymers

	return monomerCounts, polymersCreated
