from __future__ import absolute_import, division, print_function

import copy

from vivarium.core.process import Deriver
from vivarium.library.dict_utils import deep_merge
from vivarium.library.units import units
from vivarium_cell.processes.derive_globals import (
	length_from_volume,
	surface_area_from_length,
)


class WcEcoliDeriveShape(Deriver):

	defaults = {
		'width': 1,  # um
	}
	name = "wcEcoliDeriveShape"

	def __init__(self, initial_parameters=None):
		# type: (dict) -> None
		'''Derives cell length and surface area from width and volume.

		Ports:

		* **global**: Should be given the agent's boundary store.

		Arguments:
			initial_parameters (dict): A dictionary that can contain the
				follwing configuration options:

				* **width** (:py:class:`float`): Width of the cell in
				  microns
		'''
		if initial_parameters is None:
			initial_parameters = {}
		parameters = copy.deepcopy(self.defaults)
		deep_merge(parameters, initial_parameters)

		super(WcEcoliDeriveShape, self).__init__(parameters)

	def ports_schema(self):
		default_state = {
			'global': {
				'volume': 0.0 * units.fL,
				'width': self.parameters['width'],
				'length': 0.0,
				'surface_area': 0.0,
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
