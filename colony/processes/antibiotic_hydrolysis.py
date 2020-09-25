from __future__ import absolute_import, division, print_function

import copy

from vivarium.library.units import units
from vivarium.library.dict_utils import deep_merge, deep_merge_check
from vivarium_cell.processes.convenience_kinetics import ConvenienceKinetics


class AntibioticHydrolysis(ConvenienceKinetics):

	name = 'antibiotic_hydrolysis'
	defaults = {
		'kcat': 1 / units.sec,
		'Km': 1e-3 * units.mmol / units.L,
		'target': 'antibiotic',
		'initial_target_internal': 0,
		'initial_hydrolyzed_internal': 0,
		'catalyst': 'catalyst',
		'initial_catalyst': 1,
	}

	def __init__(self, initial_parameters=None):
		if initial_parameters is None:
			initial_parameters = {}
		super_defaults = super(AntibioticHydrolysis, self).defaults
		deep_merge_check(self.defaults, super_defaults)
		self.defaults.update(super_defaults)
		parameters = copy.deepcopy(self.defaults)
		deep_merge(parameters, initial_parameters)

		kcat = parameters['kcat'].to(1 / units.sec).magnitude
		km = parameters['Km'].to(units.mmol / units.L).magnitude
		hydrolyzed_key = '{}_hydrolyzed'.format(parameters['target'])

		kinetics_parameters = {
			'reactions': {
				'hydrolysis': {
					'stoichiometry': {
						('internal', parameters['target']): -1,
						('internal', hydrolyzed_key): 1,
					},
					'is_reversible': False,
					'catalyzed by': [
						('catalyst_port', parameters['catalyst'])],
				},
			},
			'kinetic_parameters': {
				'hydrolysis': {
					('catalyst_port', parameters['catalyst']): {
						'kcat_f': kcat,
						('internal', parameters['target']): km,
					},
				},
			},
			'initial_state': {
				'fluxes': {
					'hydrolysis': 0.0,
				},
				'internal': {
					parameters['target']: parameters[
						'initial_target_internal'],
					hydrolyzed_key: parameters[
						'initial_hydrolyzed_internal'],
				},
				'catalyst_port': {
					parameters['catalyst']: parameters[
						'initial_catalyst'],
				},
			},
			'port_ids': ['internal', 'catalyst_port'],
		}

		super(AntibioticHydrolysis, self).__init__(kinetics_parameters)
