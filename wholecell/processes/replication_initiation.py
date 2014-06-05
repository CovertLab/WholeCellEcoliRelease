#!/usr/bin/env python

"""
ReplicationInitiation

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process


class ReplicationInitiation(wholecell.processes.process.Process):
	""" ReplicationInitiation """

	_name = "ReplicationInitiation"

	# Constructor
	def __init__(self):


		super(ReplicationInitiation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ReplicationInitiation, self).initialize(sim, kb)

		pass

	def calculateRequest(self):
		pass

	# Calculate temporal evolution
	def evolveState(self):
		pass