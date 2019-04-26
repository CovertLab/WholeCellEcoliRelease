from __future__ import absolute_import, division, print_function

import time
import numpy as np
import math

from agent.inner import CellSimulation


TUMBLE_JITTER = 2.0 # (radians)
DEFAULT_COLOR = [color/255 for color in [76, 0 , 153]]

def add_lists(list1, list2):
	list3 = [x1 + x2 for (x1, x2) in zip(list1, list2)]
	return list3

def subtract_lists(list1, list2):
	list3 = [x1 - x2 for (x1, x2) in zip(list1, list2)]
	return list3

class Endocrine(CellSimulation):
	''' Endocrine Surrogate '''

	def __init__(self, state):
		self.initial_time = state.get('time', 0.0)
		self.local_time = 0.0
		# self.timestep = 1.0
		self.environment_change = {}
		self.volume = 1.0
		self.division_time = 100
		self.color = DEFAULT_COLOR

		# initial state
		prob_sender = 0.5
		self.sender_color = [76, 0 , 153]
		self.receiver_color_dark = [76, 0 , 153] #[0, 153 , 76]
		self.receiver_color_bright = [153, 255, 204]
		self.receiver_color_distance = subtract_lists(self.receiver_color_bright, self.receiver_color_dark)

		if np.random.uniform(0,1) < prob_sender:
			self.state = 'sender'
			self.color = [color / 255 for color in self.sender_color]
		else:
			self.state = 'receiver'
			self.color = [color/255 for color in self.receiver_color_dark]

		# self.motile_force = [0.0, 0.0] # initial magnitude and relative orientation
		self.division = []


	def update_state(self):

		if self.state == 'sender':
			self.environment_change['signal'] = 1e8
		if self.state == 'receiver':
			activation = self.external_concentrations['signal']
			if (math.exp(10000 * activation) - 1) < 1:
				brightness = [color * (math.exp(10000 * activation) - 1) for color in self.receiver_color_distance]
				activation_color = add_lists(self.receiver_color_dark, brightness)
			else:
				activation_color = self.receiver_color_bright

			self.color = [color / 255 for color in activation_color]

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
		# self.check_division()
		self.local_time = run_until

		time.sleep(1)  # pause for better coordination with Lens visualization. TODO: remove this

	def generate_inner_update(self):
		return {
			'volume': self.volume,
			# 'motile_force': self.motile_force,
			'environment_change': self.environment_change,
			'division': self.division,
			'color': self.color,
			}
