"""
Classes used to execute arbitrary code during critical parts of the simulation.
These are primarily used to support simulations that require behavior that 
cannot be modeled at this time.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/8/2014
"""

from __future__ import division

import numpy as np

class SimulationHook(object):
	_name = None

	def __init__(self):
		pass


	@classmethod
	def name(cls):
		return cls._name


	def initialize(self, sim, kb):
		pass


	def postCalcInitialConditions(self, sim):
		pass


	def preEvolveState(self, sim):
		pass


	def postEvolveState(self, sim):
		pass


	def finalize(self, sim):
		pass


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
