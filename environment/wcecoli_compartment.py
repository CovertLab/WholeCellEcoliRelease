from __future__ import absolute_import, division, print_function

import copy

from vivarium.core.experiment import Compartment

from environment.wcecoli_process import wcEcoliAgent
from environment.wcecoli_meta_division import WcEcoliMetaDivision
from environment.wcecoli_derive_shape import WcEcoliDeriveShape


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
		derive_shape = WcEcoliDeriveShape()
		return {
			'wcecoli': wcecoli_process,
			'meta_division': meta_division,
			'derive_shape': derive_shape,
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
			'derive_shape': {
				'global': boundary_path,
			},
		}
