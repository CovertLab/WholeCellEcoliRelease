from __future__ import absolute_import, division, print_function

import time
import numpy as np

from agent.inner import CellSimulation


TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [76, 0 , 153]]

class Chemotaxis(CellSimulation):
	'''
	Simple chemotaxis surrogate that can move up glucose gradients. It can take on two states 'run' and 'tumble'.
	State is a function of the current glucose concentrations, and internal CheY concentrations -- the 'memory' of
	glucose concentrations from the previous time step.

	TODO (Eran) implement mechanistic model of chemotaxis pathway. The following paper's model shoud suffice:
	Bray, Dennis, Robert B. Bourret, and Melvin I. Simon. "Computer simulation of the phosphorylation cascade controlling bacterial chemotaxis." Molecular Biology of the Cell (1993)
	'''

	def __init__(self, state):
		self.initial_time = state.get('time', 0.0)
		self.local_time = 0.0
		# self.timestep = 1.0

		self.environment_change = {}
		self.volume = 1.0
		self.division_time = 100

		# initial state
		self.state = ['tumble']
		self.external_concentrations = {
			'GLC': 0.0
		}
		self.internal_concentrations = {
			'Sensor': 5.0,  # uM, (Bray, Bourret, Simon 1993)
			'CheR': 1.0,    # uM, (Bray, Bourret, Simon 1993)
			'CheB': 2.0,    # uM, (Bray, Bourret, Simon 1993)
			'CheW': 5.0,    # uM, (Bray, Bourret, Simon 1993)
			'CheA': 5.0,    # uM, (Bray, Bourret, Simon 1993)
			'CheY': 10.0,   # uM, (Bray, Bourret, Simon 1993)
			'CheY-P': 0.0,  # uM, (Bray, Bourret, Simon 1993)
			'CheZ' : 20.0,  # uM, (Bray, Bourret, Simon 1993)
			'Motor': 0.01   # uM, (Bray, Bourret, Simon 1993)
		}
		self.motile_force = [0.0, 0.0] # initial magnitude and relative orientation
		self.division = []


	def update_state(self):
		# update state based on internal and external concentrations
		if self.external_concentrations['GLC'] >= self.internal_concentrations['CheY-P']:
			self.state = 'run'
		else:
			self.state = 'tumble'

		# update intracellular concentrations
		self.internal_concentrations['CheY-P'] = self.external_concentrations['GLC']

	def update_behavior(self):
		# update behavior based on the current state of the system
		if self.state is 'run':
			force = 0.2
			torque = 0.0
			self.motile_force = [force, torque]
		elif self.state is 'tumble':
			force = 0.05
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

		time.sleep(0.2)  # pause for better coordination with Lens visualization. TODO: remove this

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'color': DEFAULT_COLOR,
			}
