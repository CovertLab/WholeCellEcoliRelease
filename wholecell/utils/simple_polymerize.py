#!/usr/bin/env python

"""
simple_polymerize.py

A generalization of the algorithm used in our simplified (bulk) transcription
and translation submodels.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/22/2014
"""

from __future__ import division

import numpy as np

POLYMER_CREATION_FALLOFF_RATE = 0.5

def simplePolymerize(templateMonomerCounts, enzymaticLimitation,
		monomerCounts, synthesisProbabilities, randomState):

	assert np.abs(synthesisProbabilities.sum() - 1) < 1e-5

	templateLengths = templateMonomerCounts.sum(axis = 1)

	avgTemplateLength = templateLengths.mean()

	nPolymersToCreate = int(
		0.5 * min(enzymaticLimitation, monomerCounts.sum()) / avgTemplateLength
		)

	polymersCreated = np.zeros(templateLengths.size, np.int)

	nNewPolymers = np.zeros(templateLengths.size, np.int)

	canPolymerize = (synthesisProbabilities > 0)

	maxTemplateMonomerCounts = np.max(templateMonomerCounts, axis = 0)
	maxTemplateLength = np.max(templateLengths)

	while True:
		if nPolymersToCreate <= 0:
			break

		# Can only polymerize if 1) there are enough monomers and 2) there is enough enzymatic capacity
		enoughMonomers = (
			True if (monomerCounts >= maxTemplateMonomerCounts).all() else
			(monomerCounts >= templateMonomerCounts[canPolymerize, :]).all(axis = 1)
			)

		enoughEnzymaticCapacity = (
			True if (enzymaticLimitation >= maxTemplateLength) else
			(enzymaticLimitation >= templateLengths[canPolymerize])
			)

		canPolymerize[canPolymerize] = enoughMonomers & enoughEnzymaticCapacity

		if not np.any(canPolymerize):
			break

		adjustedSynthProb = canPolymerize * synthesisProbabilities
		adjustedSynthProb /= adjustedSynthProb.sum()

		# Choose numbers of each polymer to synthesize
		nNewPolymers[:] = 0
		nNewPolymers[canPolymerize] = randomState.multinomial(
			nPolymersToCreate,
			adjustedSynthProb[canPolymerize]
			)

		# TODO: determine why this assertion fails on Hudson
		# assert nNewPolymers.sum() == nPolymersToCreate

		nonzeroNewPolymers = (nNewPolymers > 0)
		monomersUsed = np.dot(
			nNewPolymers[nonzeroNewPolymers],
			templateMonomerCounts[nonzeroNewPolymers, :]
			)

		enzymaticCapacityUsed = monomersUsed.sum()

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

		nPolymersToCreate = int(nPolymersToCreate / POLYMER_CREATION_FALLOFF_RATE)

	return monomerCounts, polymersCreated
