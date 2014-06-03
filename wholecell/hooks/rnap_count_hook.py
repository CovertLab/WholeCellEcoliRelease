"""
Hook that fixes the RNA poly level to grow exponentially.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import division

import numpy as np

from wholecell.hooks.hook import SimulationHook

class RnapCountHook(SimulationHook):
	_name = "RnapCountHook"

	_rnapSubunitIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]
	_rnapSubunitStoich = np.array([2, 1, 1, 1])

	_activeRnapName = "activeRnaPoly"

	_verbose = False

	def initialize(self, sim, kb):
		self._cellCycleLen = kb.cellCycleLen.to('s').magnitude

		self._bulkMolecules = sim.states["BulkMolecules"]
		self._uniqueMolecules = sim.states["UniqueMolecules"]

		self._initialFreeRnap = None


	def postCalcInitialConditions(self, sim):
		self._initialFreeRnap = self._nTotalRnap()


	def preEvolveState(self, sim):
		nTotalRnap = self._nTotalRnap()

		nExpectedRnap = np.int(self._initialFreeRnap *
			np.exp(np.log(2)/self._cellCycleLen * sim.time()))

		nRnapDifference = nExpectedRnap - nTotalRnap

		freeRnapSubunits = self._bulkMolecules.container.countsView(self._rnapSubunitIds)

		if nRnapDifference > 0:
			if self._verbose:
				print "Increased RNAP subunit availability by {}".format(nRnapDifference)

			freeRnapSubunits.countsInc(
				nRnapDifference * self._rnapSubunitStoich
				)

		elif nRnapDifference < 0:
			nToRemove = min(
				-nRnapDifference,
				freeRnapSubunits.counts().min()
				)

			if self._verbose:
				print "Decreased RNAP subunit availability by {}".format(-nRnapDifference)

			freeRnapSubunits.countsDec(
				nToRemove * self._rnapSubunitStoich
				)


	def _nTotalRnap(self):
		freeRnapSubunits = self._bulkMolecules.container.counts(self._rnapSubunitIds)

		nMaxFreeRnap = np.min(freeRnapSubunits // self._rnapSubunitStoich)

		activeRnap = self._uniqueMolecules.container.objectsInCollection(self._activeRnapName)

		nActiveRnap = len(activeRnap)

		return nMaxFreeRnap + nActiveRnap
