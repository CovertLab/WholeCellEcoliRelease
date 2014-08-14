#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process

class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	# Constructor
	def __init__(self):
		super(PolypeptideElongation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(PolypeptideElongation, self).initialize(sim, kb)

		# Load constants

		# Create views on state


	def calculateRequest(self):
		pass


	def evolveState(self):
		pass
