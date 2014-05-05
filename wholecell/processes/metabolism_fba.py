#!/usr/bin/env python

"""
MetabolismFba

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/5/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class MetabolismFba(wholecell.processes.process.Process):
	""" MetabolismFba """

	_name = "MetabolismFba"

	# Constructor
	def __init__(self):
		self.lowerBounds = None
		self.upperBounds = None
		self.stoichMatrix = None
		self.objective = None

		super(MetabolismFba, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(MetabolismFba, self).initialize(sim, kb)
		print 'init'


	def calculateRequest(self):
		print 'request'


	# Calculate temporal evolution
	def evolveState(self):
		print 'evolve'
