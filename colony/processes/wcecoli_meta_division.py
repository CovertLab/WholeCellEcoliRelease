from __future__ import absolute_import, division, print_function

import uuid
import copy

from vivarium.core.process import Deriver
from vivarium.library.dict_utils import deep_merge


def divide_empty_list(state):
	'''Assign each daughter variable an empty list
	'''
	return [[], []]


class WcEcoliMetaDivision(Deriver):

	defaults = {
		'daughter_path': tuple(),
		'id_function': lambda: str(uuid.uuid1()),
	}
	name = 'wcEcoliMetaDivision'

	def __init__(self, initial_parameters):
		'''Implement division for wcEcoli simulations

        When wcEcoli triggers division internally, it passes through its
        update a list of 2 paths to states for the daughter cells. The
        :py:class:`colony.processes.wcecoli.WcEcoli` process stores this
        in the ``division`` variable of the ``global`` port, which is
        associated with the same store as the ``global`` port for this
        process. When this deriver runs then, it uses these paths to
        instantiate new WcEcoli processes and uses the ``_divide``
        operation to trigger the division in the hierarchy.

		Ports:

		* ``global``: The store with the ``division`` variable.
		* ``cells``: The store that contains the roots of each agent's
		  subtree in the hierarchy. Typically, the store is called
		  ``agents``.

		Arguments:
			initial_parameters(dict): Accepts the following
				configuration keys:

				* **id_function**: Function that generates a unique
				  identifier for the agents generated upon
				  division.
				* **compartment**
				  (:py:class:`vivarium.core.experiment.Compartment`):
				  The compartment object that will generate the agents.
				* **daughter_path** (:py:class:`str`): Hierarchy path to
				  the store where the agents will be created, relative
				  to the store associated with the``cells`` port.
				* **agent_id** (:py:class:`str`): The identifier for the
				  mother agent.

				Unlike many other processes, this configuration
				dictionary MUST be provided and MUST specify
				``compartment`` and ``agent_id``.
		'''
		# must provide a compartment to generate new daughters
		assert 'compartment' in initial_parameters
		# must provide an agent ID
		assert 'agent_id' in initial_parameters
		parameters = copy.deepcopy(self.defaults)
		deep_merge(parameters, initial_parameters)

		self.id_function = parameters['id_function']
		self.compartment = parameters['compartment']
		self.daughter_path = parameters['daughter_path']
		self.agent_id = parameters['agent_id']

		super(WcEcoliMetaDivision, self).__init__(parameters)

	def ports_schema(self):
		return {
			'global': {
				'division': {
					'_default': [],
					'_updater': 'set',
					'_divider': divide_empty_list,
				},
			},
			'cells': {
				'*': {},
			}
		}

	def next_update(self, timestep, states):
		daughter_configs = states['global']['division']
		if daughter_configs == []:
			return {}

		daughter_updates = []

		for daughter_config in daughter_configs:
			daughter_id = daughter_config['id']
			assert daughter_id != self.agent_id
			compartment = self.compartment.generate({
				'agent_id': daughter_id,
				'start_time': daughter_config['start_time'],
				'volume': daughter_config['volume'],
				'agent_config': {
					'files': [daughter_config['inherited_state_path']],
				},
				'seed': daughter_config['seed'],
			})
			daughter_updates.append({
				'daughter': daughter_id,
				'path': (daughter_id,) + self.daughter_path,
				'processes': compartment['processes'],
				'topology': compartment['topology'],
				'initial_state': {},
			})

		return {
			'cells': {
				'_divide': {
					'mother': self.agent_id,
					'daughters': daughter_updates,
				},
			},
		}
