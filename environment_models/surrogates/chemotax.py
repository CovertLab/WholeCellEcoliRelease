from __future__ import absolute_import, division, print_function

import time
import numpy as np
import random

from agent.inner import CellSimulation

class Chemotax(CellSimulation):


	def __init__(self):
		self.local_time = 0
		self.concentrations = {}

		self.volume = 1.0
		self.division = False
		self.environment_change = {}


	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		self.concentrations = update['concentrations']

		self.environment_change = {}
		for molecule in self.concentrations.iterkeys():
			self.environment_change[molecule] = 0

	def run_incremental(self, run_until):
		self.local_time = run_until

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			# 'division': self.division,
			'environment_change': self.environment_change,
			}

	def finalize(self):
		pass
