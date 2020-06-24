from __future__ import absolute_import, division, print_function

import copy

from vivarium.core.process import Deriver
from vivarium.library.units import units
from vivarium.library.dict_utils import deep_merge
from vivarium.processes.derive_globals import (
	length_from_volume,
	surface_area_from_length,
)


class WcEcoliDeriveGlobals(Deriver):
	"""
	Process for deriving volume, mmol_to_counts, and shape from the cell mass
	"""

	defaults = {
		'width': 1,  # um
		'volume_units': units.fL,
		'length_units': units.um,  # Same as width units
	}

	def __init__(self, initial_parameters=None):
		if initial_parameters is None:
			initial_parameters = {}
		parameters = copy.deepcopy(self.defaults)
		deep_merge(parameters, initial_parameters)

		self.volume_units = parameters['volume_units']
		self.length_units = parameters['length_units']
		self.width = parameters['width'] * self.length_units

		ports = {
			'global': [
				'volume',
				'width',
				'length',
				'surface_area',
			],
		}

		parameters = {}
		parameters.update(initial_parameters)

		super(WcEcoliDeriveGlobals, self).__init__(ports, parameters)

	def ports_schema(self):
		default_state = {
			'global': {
				'volume': 0.0 * self.volume_units,
				'width': self.width,
				'length': 0.0 * self.length_units,
				'surface_area': 0.0 * self.length_units ** 2,
			}
		}

		schema = {
			'global': {
				variable: {
					'_updater': 'set',
					'_emit': True,
					'_divider': (
						'set' if variable == 'width' else 'split'
					),
					'_default': default_state['global'][variable]
				}
				for variable in default_state['global']
			}
		}
		return schema

	def next_update(self, timestep, states):
		width = states['global']['width']
		volume = states['global']['volume']

		length = length_from_volume(volume.magnitude, width)
		surface_area = surface_area_from_length(length, width)

		return {
			'global': {
				'length': length,
				'surface_area': surface_area,
			},
		}
