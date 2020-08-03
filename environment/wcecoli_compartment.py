from __future__ import absolute_import, division, print_function

import copy

from vivarium.core.process import Generator
from vivarium.processes.death import DeathFreezeState
from vivarium.processes.antibiotic_transport import AntibioticTransport
from vivarium.processes.diffusion_cell_environment import (
	CellEnvironmentDiffusion,
)
from vivarium.processes.derive_concentrations import (
	DeriveConcentrations)

from environment.wcecoli_process import wcEcoliAgent
from environment.wcecoli_meta_division import WcEcoliMetaDivision
from environment.wcecoli_derive_shape import WcEcoliDeriveShape


INITIAL_INTERNAL_ANTIBIOTIC = 0
INITIAL_EXTERNAL_ANTIBIOTIC = 1
ANTIBIOTIC_KEY = 'rifampicin'
PUMP_KEY = 'TRANS-CPLX-201[s]'
PORIN_KEY = 'porin'


class WcEcoliCell(Generator):

	defaults = {
		'boundary_path': ('boundary',),
		'agents_path': ('..', '..', 'agents'),
		'fields_path': ('..', '..', 'fields'),
		'dimensions_path': ('..', '..', 'dimensions'),
		'antibiotic_transport': {
			'initial_pump': 0.0,
			'initial_internal_antibiotic': INITIAL_INTERNAL_ANTIBIOTIC,
			'intial_external_antibiotic': INITIAL_EXTERNAL_ANTIBIOTIC,
			'pump_kcat': 1,
			'pump_key': PUMP_KEY,
			'antibiotic_key': ANTIBIOTIC_KEY,
		},
		'death': {
			'detectors': {
				'antibiotic': {
					'antibiotic_threshold': 0.86,
					'antibiotic_key': ANTIBIOTIC_KEY,
				},
			},
			'targets': [
				'wcecoli', 'meta_division', 'death',
				'antibiotic_transport'
			],
		},
		'derive_shape': {},
		'cell_environment_diffusion': {
			'default_state': {
				'membrane': {
					PORIN_KEY: 1e-7,
				},
				'external': {
					ANTIBIOTIC_KEY: INITIAL_EXTERNAL_ANTIBIOTIC,
				},
				'internal': {
					ANTIBIOTIC_KEY: INITIAL_INTERNAL_ANTIBIOTIC,
				},
			},
			'molecules_to_diffuse': [ANTIBIOTIC_KEY],
			'permeabilities': {
				PORIN_KEY: 1e-11,
			},
		},
		'derive_colony_shape': {
			'alpha': 0.5,
		},
		'derive_concentrations': {
			'concentration_keys': [PUMP_KEY],
		},
        '_schema': {
            'cell_environment_diffusion': {
                'membrane': {
                    PORIN_KEY: {
                        '_divider': 'set',
                    },
                },
            },
        },
	}

	def generate_processes(self, config):
		wcecoli_process = wcEcoliAgent(config)
		meta_division = WcEcoliMetaDivision({
			'agent_id': config['agent_id'],
			'compartment': self,
		})
		derive_shape = WcEcoliDeriveShape(config['derive_shape'])
		death = DeathFreezeState(config['death'])
		antibiotic_transport = AntibioticTransport(
			config['antibiotic_transport'])
		cell_environment_diffusion = CellEnvironmentDiffusion(
			config['cell_environment_diffusion'])
		derive_concentrations = DeriveConcentrations(
			config['derive_concentrations'])
		return {
			'derive_concentrations': derive_concentrations,
			'wcecoli': wcecoli_process,
			'meta_division': meta_division,
			'derive_shape': derive_shape,
			'death': death,
			'antibiotic_transport': antibiotic_transport,
			'cell_environment_diffusion': cell_environment_diffusion,
		}

	def generate_topology(self, config=None):
		boundary_path = config['boundary_path']
		return {
			'derive_concentrations': {
				'global': boundary_path,
				'counts': boundary_path + ('bulk_molecules_report',),
				'concentrations': (
					boundary_path + ('bulk_molecule_concentrations',)
				),
			},
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
				'dimensions': config['dimensions_path'],
			},
			'meta_division': {
				'global': boundary_path,
				'cells': config['agents_path'],
			},
			'derive_shape': {
				'global': boundary_path,
			},
			'death': {
				'global': boundary_path,
				'internal': boundary_path + ('cytoplasm',),
				'processes': tuple(),
			},
			'antibiotic_transport': {
				'internal': boundary_path + ('cytoplasm',),
				'external': boundary_path + ('external',),
				'pump_port': (
					boundary_path + ('bulk_molecule_concentrations',)
				),
				'fields': config['fields_path'],
				'fluxes': boundary_path + ('fluxes',),
				'global': boundary_path,
				'dimensions': config['dimensions_path'],
			},
			'cell_environment_diffusion': {
				'membrane': boundary_path + ('membrane',),
				'internal': boundary_path + ('cytoplasm',),
				'external': boundary_path + ('external',),
				'fields': config['fields_path'],
				'global': boundary_path,
				'dimensions': config['dimensions_path'],
			},
		}
