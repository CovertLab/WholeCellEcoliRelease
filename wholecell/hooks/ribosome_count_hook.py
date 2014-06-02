"""
Hook that fixes the ribosome level to grow exponentially.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import division

import numpy as np

from wholecell.hooks.hook import SimulationHook

class RibosomeCountHook(SimulationHook):
	_name = "RibosomeCountHook"

	_ribosomeSubunitIds = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]
	_ribosomeSubunitStoich = np.array([1, 1, 1])

	_activeRibosomeId = "activeRibosome"

	_verbose = False

	def initialize(self, sim, kb):
		self._cellCycleLen = kb.cellCycleLen.to('s').magnitude

		self._bulkMolecules = sim.states["BulkMolecules"]
		self._uniqueMolecules = sim.states["UniqueMolecules"]

		self._initialFreeRibosomes = None


	def postCalcInitialConditions(self, sim):
		self._initialFreeRibosomes = self._nTotalRibosome()


	def preEvolveState(self, sim):
		nTotalRibosome = self._nTotalRibosome()

		nExpectedRibosome = np.int(self._initialFreeRibosomes *
			np.exp(np.log(2)/self._cellCycleLen * sim.time()))

		nRibDifference = nExpectedRibosome - nTotalRibosome

		freeRibSubunits = self._bulkMolecules.container.countsView(self._ribosomeSubunitIds)

		if nRibDifference > 0:
			if self._verbose:
				print "Increased ribosome subunit availability by {}".format(nRibDifference)

			freeRibSubunits.countsInc(
				nRibDifference * self._ribosomeSubunitStoich
				)

		elif nRibDifference < 0:
			nToRemove = min(
				-nRibDifference,
				freeRibSubunits.counts().min()
				)

			if self._verbose:
				print "Decreased ribosome subunit availability by {}".format(-nRibDifference)

			freeRibSubunits.countsDec(
				nToRemove * self._ribosomeSubunitStoich
				)


	def _nTotalRibosome(self):
		freeRibSubunits = self._bulkMolecules.container.counts(self._ribosomeSubunitIds)

		nMaxFreeRibs = np.min(freeRibSubunits // self._ribosomeSubunitStoich)

		activeRib = self._uniqueMolecules.container.objectsInCollection(self._activeRibosomeId)

		nActiveRib = len(activeRib)

		return nMaxFreeRibs + nActiveRib
