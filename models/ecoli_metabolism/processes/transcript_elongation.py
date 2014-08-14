#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process

class TranscriptElongation(wholecell.processes.process.Process):
	""" TranscriptElongation """

	_name = "TranscriptElongation"

	# Constructor
	def __init__(self):
		super(TranscriptElongation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(TranscriptElongation, self).initialize(sim, kb)

		# Load constants

		# Create views on state


	def calculateRequest(self):
		pass


	def evolveState(self):
		pass
