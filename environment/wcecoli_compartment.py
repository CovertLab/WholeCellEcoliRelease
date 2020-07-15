from __future__ import absolute_import, division, print_function

import copy

from vivarium.core.process import Generator
from vivarium.processes.death import DeathFreezeState
from vivarium.processes.injector import Injector
from vivarium.processes.exchange_suppressor import ExchangeSuppressor

from environment.wcecoli_process import wcEcoliAgent
from environment.wcecoli_meta_division import WcEcoliMetaDivision
from environment.wcecoli_derive_shape import WcEcoliDeriveShape


class WcEcoliCell(Generator):

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
		death_config = {
			'detectors': {
				'antibiotic': {
					'antibiotic_threshold': 10.0,
					'antibiotic_key': 'rifampicin',
				},
			},
			'targets': ['wcecoli', 'meta_division', 'death', 'injector'],
		}
		death = DeathFreezeState(death_config)
		injector_config = {
			'substrate_rate_map': {
				'rifampicin': 1.0,
			},
		}
		injector = Injector(injector_config)
		exchange_suppressor = ExchangeSuppressor()
		return {
			'wcecoli': wcecoli_process,
			'meta_division': meta_division,
			'derive_shape': derive_shape,
			'death': death,
			'injector': injector,
			'exchange_suppressor': exchange_suppressor,
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
			'death': {
				'global': boundary_path,
				'internal': boundary_path + ('cytoplasm',),
			},
			'injector': {
				'internal': boundary_path + ('cytoplasm',),
			},
			'exchange_suppressor': {
				'exchange': boundary_path + ('exchange',),
				'global': boundary_path,
			},
		}
