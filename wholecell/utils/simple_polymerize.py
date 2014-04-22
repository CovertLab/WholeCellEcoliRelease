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

OVERALL_SYNTHESIS_PROBABILITY_MINIMUM = 1e-3

def simplePolymerize(templateMonomerCounts, enzymaticLimitation,
		monomerCounts, synthesisProbabilities, randStream):

	templateLengths = templateMonomerCounts.sum(1)

	avgTemplateLength = templateLengths.mean()

	nPolymersToCreate = int(enzymaticLimitation / avgTemplateLength)
	
	polymersCreated = np.zeros_like(templateLengths)

	while enzymaticLimitation > 0 and nPolymersToCreate > 0:
		# Can only polymerize if 1) there are enough monomers and 2) there is enough enzymatic power
		canPolymerize = ((monomerCounts >= templateMonomerCounts).all(1) &
			(enzymaticLimitation >= templateLengths))

		if not np.any(canPolymerize):
			break

		# Break if overall synthesis probability is too low
		if (synthesisProbabilities[canPolymerize].sum() <
				OVERALL_SYNTHESIS_PROBABILITY_MINIMUM):
			break

		# TODO: Give unpolymerizable templates a synthesis probability of zero

		# Choose numbers of each polymer to synthesize
		nNewPolymers = randStream.mnrnd(nPolymersToCreate, synthesisProbabilities)

		monomersUsedByTemplate = np.dot(nNewPolymers, templateMonomerCounts)
		monomersUsed = monomersUsedByTemplate.sum(0)

		enzymaticPowerUsed = np.dot(nNewPolymers, templateLengths).sum()

		# Reduce number of polymers to create if not enough monomers
		if (monomerCounts < monomersUsed).any():
			if nPolymersToCreate > 0:
				nPolymersToCreate //= 2
				continue

			else:
				break

		# Reduce number of polymers to create if not enough enzymatic power
		if enzymaticLimitation < enzymaticPowerUsed:
			if nPolymersToCreate > 0:
				nPolymersToCreate //= 2
				continue

			else:
				break

		# Use resources and create the polymers

		enzymaticLimitation -= enzymaticPowerUsed

		monomerCounts -= monomersUsed

		polymersCreated += nNewPolymers

	return monomerCounts, polymersCreated
