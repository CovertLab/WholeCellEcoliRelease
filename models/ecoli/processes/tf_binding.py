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

		self.pTfBound = sim_data.process.transcription_regulation.pTfBound
		self.tfKd = sim_data.process.transcription_regulation.tfKd
		self.tfNTargets = sim_data.process.transcription_regulation.tfNTargets

		self.alphaView = self.bulkMoleculesView(alphaNames)
		self.tfBoundViews = {}
		self.tfMoleculeViews = {}
		for tf in self.tfs:
			self.tfBoundViews[tf] = self.bulkMoleculesView(tfNames[tf])
			self.tfMoleculeViews[tf] = self.bulkMoleculeView(tf + "[c]")


	def calculateRequest(self):
		self.alphaView.requestAll()
		for view in self.tfMoleculeViews.itervalues():
			view.requestAll()
		for view in self.tfBoundViews.itervalues():
			view.requestAll()


	def evolveState(self):
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)
		cellVolume = cellMass / self.cellDensity
		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		self.alphaView.countsIs(1)

		for tf in self.tfs:
			tfFreeCounts = self.tfMoleculeViews[tf].count()
			tfBoundCounts = self.tfBoundViews[tf].counts()
			if tfFreeCounts == 0:
				continue
			tfKd = self.tfKd[tf]
			promoterConc = countsToMolar * self.tfNTargets[tf]
			tfConc = countsToMolar * tfFreeCounts

			pTfBound = self.pTfBound(
				tfKd.asNumber(units.nmol / units.L),
				promoterConc.asNumber(units.nmol / units.L),
				tfConc.asNumber(units.nmol / units.L)
				)
			nToBind = int(stochasticRound(self.randomState, tfFreeCounts * pTfBound))
			if nToBind == 0:
				continue

			self.tfMoleculeViews[tf].countDec(nToBind)

			nBound = tfBoundCounts.sum()
			boundLocs = np.zeros_like(tfBoundCounts)
			boundLocs[
				self.randomState.choice(tfBoundCounts.size, size = np.min((nToBind, tfBoundCounts.size)), replace = False)
				] = 1

			self.tfBoundViews[tf].countsIs(boundLocs)
			self.tfMoleculeViews[tf].countInc(nBound)