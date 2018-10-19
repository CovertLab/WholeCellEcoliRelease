from __future__ import absolute_import, division, print_function

import time
import numpy as np
import random

from agent.inner import CellSimulation


TUMBLE_JITTER = 0.4 # (radians)


class Chemotaxis(CellSimulation):

	def __init__(self):
		self.initial_time = 0.0
		self.local_time = 0.0
		self.timestep = 1.0
		self.external_concentrations = {'GLC[p]': 0.0}
		self.environment_change = {}
		self.volume = 1.0

		self.memory = {'GLC[p]': 0.0}
		self.state = None
		self.motile_force = [0.0, 0.0] # magnitude and relative orientation

		self.division = []
		self.division_time = 500

		# TODO (Eran) include state names (['run', 'tumble']) and make surrogate into a proper state machine


	def update_state(self):

		# update state based on change in concentration
		if self.external_concentrations['GLC[p]'] >= self.memory['GLC[p]']:
			# run
			self.motile_force = [0.02, 0.0]
		else:
			# tumble
			self.motile_force = [0.005, np.random.normal(scale=TUMBLE_JITTER)]

		# update memory to current concentrations
		self.memory = self.external_concentrations

		# divide cell
		if self.local_time >= self.initial_time + self.division_time:
			self.divide()

		time.sleep(0.1) # pause for better coordination with Lens visualization. TODO: remove this

	def divide(self):
		self.division = [{'time': self.local_time}, {'time': self.local_time}]

		return self.division

	def time(self):
		return self.local_time

	def apply_outer_update(self, update):
		self.external_concentrations = update['concentrations']

		self.environment_change = {}
		for molecule in self.external_concentrations.iterkeys():
			self.environment_change[molecule] = 0

	def run_incremental(self, run_until):
		# update state once per message exchange
		self.update_state()
		self.local_time = run_until

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			}

	def synchronize_state(self, state):
		if 'time' in state:
			self.initial_time = state['time']