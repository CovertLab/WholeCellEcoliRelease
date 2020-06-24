from __future__ import absolute_import, division, print_function

import copy

from vivarium.core.experiment import Compartment
from wcecoli_process import wcEcoliAgent
from wcecoli_meta_division import WcEcoliMetaDivision
from wcecoli_derive_globals import WcEcoliDeriveGlobals


class WcEcoliCell(Compartment):

	defaults = {
		'boundary_path': ('boundary',),
		'agents_path': ('..', '..', 'agents'),
	}

	def __init__(self, config=None):
		if config is None:
			config = {}
		self.config = copy.deepcopy(self.defaults)
		self.config.update(config)

	def generate_processes(self, config):
		wcecoli_process = wcEcoliAgent(config)
		meta_division = WcEcoliMetaDivision({
			'agent_id': config['agent_id'],
			'compartment': self,
		})
		derive_globals = WcEcoliDeriveGlobals()
		return {
			'wcecoli': wcecoli_process,
			'meta_division': meta_division,
			'derive_globals': derive_globals
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
			},
			'meta_division': {
				'global': boundary_path,
				'cells': self.config['agents_path'],
			},
			'derive_globals': {
				'global': boundary_path,
			},
		}
