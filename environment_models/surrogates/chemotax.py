from __future__ import absolute_import, division, print_function

import time
import numpy as np
import random

from agent.inner import CellSimulation


TUMBLE_JITTER = 0.4 # (radians/s)


class Chemotax(CellSimulation):

	def __init__(self):
		self.local_time = 0
		self.timestep = 1.0
		self.external_concentrations = {'GLC[p]': 0.0}
		self.environment_change = {}
		self.volume = 1.0

		self.memory = {'GLC[p]': 0.0}
		self.state = None
		self.motile_force = [0.0, 0.0] # magnitude and relative orientation

		# self.division = False
		# TODO (Eran) include state names (['run', 'tumble'])


	def update_state(self):

		# update state based on change in concentration
		if self.external_concentrations['GLC[p]'] >= self.memory['GLC[p]']:
			# run
			self.motile_force = [0.1, 0.0]
		else:
			# tumble
			self.motile_force = [0.02, np.random.normal(scale=TUMBLE_JITTER)]

		# update memory to current concentrations
		self.memory = self.external_concentrations

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		self.external_concentrations = update['concentrations']

		self.environment_change = {}
		for molecule in self.external_concentrations.iterkeys():
			self.environment_change[molecule] = 0

	def run_incremental(self, run_until):
		# update state once per update
		self.update_state()
		self.local_time = run_until

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			}
