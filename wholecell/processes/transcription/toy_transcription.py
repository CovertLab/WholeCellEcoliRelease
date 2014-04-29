#!/usr/bin/env python

"""
ToyTranscription

Creates unique RNA polymerases and binds them to an imaginary chromosome.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/28/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class ToyTranscription(wholecell.processes.process.Process):
	""" ToyTranscription """

	_name = "ToyTranscription"

	# Constructor
	def __init__(self):
		self.polymerizationRate = 48

		super(ToyTranscription, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyTranscription, self).initialize(sim, kb)

		self.rnaPolyRegions = self.chromosomeMoleculesView('RNA polymerase',
			self.polymerizationRate, 0, False)


	def calculateRequest(self):
		self.rnaPolyRegions.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		print len(self.rnaPolyRegions.regions())
