from __future__ import absolute_import, division, print_function

import copy

from vivarium.core.experiment import Compartment
from wcecoli_process import wcEcoliAgent


class WcEcoliCell(Compartment):

	defaults = {
		'boundary_path': ('boundary',),
	}

	def __init__(self, config=None):
		if config is None:
			config = {}
		self.config = self.defaults
		self.config.update(config)

	def generate_processes(self, config=None):
		if config is None:
			config = {}
		wcecoli_process = wcEcoliAgent(config)
		return {
			'wcecoli': wcecoli_process,
		}

	def generate_topology(self, config=None):
		boundary_path = self.config['boundary_path']
		return {
			'wcecoli': {
				'bulk_molecules_report': (
					boundary_path + ('bulk_molecules_report',)
				),
				'unique_molecules_report': (
					boundary_path + ('unique_molecules_report',)
				),
				'listeners_report': (
					boundary_path + ('listeners_report',)
				),
				'global': boundary_path,
				'external': (
					boundary_path + ('external',)
				),
				'exchange': boundary_path + ('exchange',),
			}
		}
