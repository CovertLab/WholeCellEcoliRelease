#!/usr/bin/env python

"""
Partition

Partition base class.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/10/2014
"""

class Partition(object):
	_state = None
	_process = None

	def __init__(self, state, process):
		self._state = state
		self._process = process


	def state(self):
		return self._state


	def process(self):
		return self._process


	def allocate(self):
		pass


	def request(self):
		raise NotImplementedError()
