from __future__ import absolute_import, division, print_function

import time
import numpy as np
import random

from agent.inner import CellSimulation


TUMBLE_JITTER = 0.4 # (radians)

class Chemotaxis(CellSimulation):
	'''
	Simple chemotaxis surrogate that can move up glucose gradients. It can take on two states 'run' and 'tumble'.
	State is a function of the current glucose concentrations, and internal CheY concentrations -- the 'memory' of
	glucose concentrations from the previous time step.

	TODO (Eran) make surrogate into a proper state machine
	TODO (Eran) make transition rates between states a function of CheY
	TODO (Eran) dynamically update CheY
	TODO (Eran) fit transition rates to experimental results
	'''

	def __init__(self):
		self.initial_time = 0.0
		self.local_time = 0.0
		# self.timestep = 1.0
		self.environment_change = {}
		self.volume = 1.0
		self.division_time = 100

		# initial state
		self.state = ['tumble']
		self.external_concentrations = {
			'GLC[p]': 0.0
		}
		self.internal_concentrations = {
			'CheY': 0.0,
			'CheY-P': 0.0,
			'CheZ' : 0.0,
			'CheA' : 0.0,
		}
		self.motile_force = [0.0, 0.0] # initial magnitude and relative orientation
		self.division = []


	def update_state(self):
		# update state based on internal and external concentrations

		if self.external_concentrations['GLC[p]'] >= self.internal_concentrations['CheY-P']:
			self.state = 'run'
		else:
			self.state = 'tumble'

		# update intracellular concentrations
		self.internal_concentrations['CheY-P'] = self.external_concentrations['GLC[p]']

	def update_behavior(self):
		# update behavior based on the current state of the system

		if self.state is 'run':
			force = 0.02
			torque = 0.0
			self.motile_force = [force, torque]
		elif self.state is 'tumble':
			force = 0.005
			torque = np.random.normal(scale=TUMBLE_JITTER)
			self.motile_force = [force, torque]

	def check_division(self):
		# update division state based on time since initialization

		if self.local_time >= self.initial_time + self.division_time:
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
		self.update_behavior()
		# self.check_division()
		self.local_time = run_until

		time.sleep(1.0)  # pause for better coordination with Lens visualization. TODO: remove this

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
