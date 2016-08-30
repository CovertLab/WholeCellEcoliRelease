#!/usr/bin/env python

"""
TfBinding

Bind transcription factors to DNA

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/14/16
"""

from __future__ import division

import warnings

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

class TfBinding(wholecell.processes.process.Process):
	""" TfBinding """

	_name = "TfBinding"

	# Constructor
	def __init__(self):
		super(TfBinding, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TfBinding, self).initialize(sim, sim_data)

		recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
		self.tfs = sorted(set([x.split("__")[-1] for x in recruitmentColNames if x.split("__")[-1] != "alpha"]))
		alphaNames = [x for x in recruitmentColNames if x.endswith("__alpha")]
		tfNames = {}
		for tf in self.tfs:
			tfNames[tf] = [x for x in recruitmentColNames if x.endswith("__" + tf)]


		self.nAvogadro = sim_data.constants.nAvogadro
		self.cellDensity = sim_data.constants.cellDensity

		self.tfNTargets = sim_data.process.transcription_regulation.tfNTargets
		self.pPromoterBoundTF = sim_data.process.transcription_regulation.pPromoterBoundTF

		self.alphaView = self.bulkMoleculesView(alphaNames)
		self.tfBoundViews = {}
		self.tfMoleculeActiveView = {}
		self.tfMoleculeInactiveView = {}
		for tf in self.tfs:
			self.tfBoundViews[tf] = self.bulkMoleculesView(tfNames[tf])
			self.tfMoleculeActiveView[tf] = self.bulkMoleculeView(tf + "[c]")
			self.tfMoleculeInactiveView[tf] = self.bulkMoleculeView(sim_data.process.equilibrium.getUnbound(tf + "[c]"))


	def calculateRequest(self):
		self.alphaView.requestAll()
		for view in self.tfMoleculeActiveView.itervalues():
			view.requestAll()
		for view in self.tfBoundViews.itervalues():
			view.requestAll()


	def evolveState(self):
		self.alphaView.countsIs(1)

		nTfs = len(self.tfs)
		pPromotersBound = np.zeros(nTfs, np.float64)
		nPromotersBound = np.zeros(nTfs, np.float64)
		nActualBound = np.zeros(nTfs, np.float64)

		for i, tf in enumerate(self.tfs):
			tfActiveCounts = self.tfMoleculeActiveView[tf].count()
			tfInactiveCounts = self.tfMoleculeInactiveView[tf].total()
			tfBoundCounts = self.tfBoundViews[tf].counts()
			promoterCounts = self.tfNTargets[tf]

			self.tfBoundViews[tf].countsIs(0)
			self.tfMoleculeActiveView[tf].countInc(tfBoundCounts.sum())

			if tfActiveCounts == 0:
				continue

			pPromoterBound = self.pPromoterBoundTF(tfActiveCounts, tfInactiveCounts)
			nToBind = int(stochasticRound(self.randomState, promoterCounts * pPromoterBound))

			if nToBind == 0:
				continue

			boundLocs = np.zeros_like(tfBoundCounts)
			boundLocs[
				self.randomState.choice(tfBoundCounts.size, size = np.min((nToBind, tfBoundCounts.size, self.tfMoleculeActiveView[tf].count())), replace = False)
				] = 1

			self.tfMoleculeActiveView[tf].countDec(boundLocs.sum())
			self.tfBoundViews[tf].countsIs(boundLocs)

			pPromotersBound[i] = pPromoterBound
			nPromotersBound[i] = pPromoterBound * self.tfNTargets[tf]
			nActualBound[i] = boundLocs.sum()

		self.writeToListener("RnaSynthProb", "pPromoterBound", pPromotersBound)
		self.writeToListener("RnaSynthProb", "nPromoterBound", nPromotersBound)
		self.writeToListener("RnaSynthProb", "nActualBound", nActualBound)
