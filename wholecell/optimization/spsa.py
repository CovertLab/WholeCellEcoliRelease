"""
Simultaneous Perturbation Stochastic Approximation (SPSA) algorithm.

Randomly perturbs all parameters up or down to calculate a gradient from two
simulations.
"""

from typing import Dict, List, Tuple

import numpy as np

from wholecell.optimization.base_solver import BaseSolver


class SPSA(BaseSolver):
	def __init__(self, method, args):
		super().__init__(method, args)

		self.alpha = args.alpha
		self.gamma = args.gamma

	### Inherited method implementations ###

	def get_parameter_updates(self, original_values: Dict, objectives: List[float], paths: List[str]) -> Dict:
		"""
		Get the new parameter values based on the results of the current
		iteration.

		Args:
			original_values: the original values of parameters at the start of
				the iteration
			objectives: objective values for each variant run in this iteration
			paths: paths to the data objects (raw_data or sim_data) for each
				variant in this iteration

		Returns:
			mapping of each parameter to be modified to its update (the amount
			to change the parameter by)
		"""

		updates = {}
		at, _ = self.get_spsa_params()

		for param, original_value in original_values.items():
			objective_diff = objectives[1] - objectives[0]
			parameter_diff = self.get_param(param, paths[1]) - self.get_param(param, paths[0])

			# Get similar scaling for each parameter even if they span a wide range of magnitudes
			# and gets the update in the correct units if the parameter has units
			parameter_scaling = original_value**2
			updates[param] = -at * objective_diff / parameter_diff * parameter_scaling

		return updates

	def get_parameter_perturbations(self, index: int) -> Tuple[Dict, Dict]:
		"""
		New parameter values for the given parameter set index within this
		iteration.

		Args:
			index: relative variant within the current iteration (0 is the first
				paramter set, 1 is the second and so on) so that parameters can
				be modified in different ways for each parameter set

		Returns:
			mapping of each parameter to be modified to its new value
		"""

		raw_data_perturbations = {}
		sim_data_perturbations = {}

		if index == 0:
			direction = 1
		elif index == 1:
			direction = -1
		else:
			raise ValueError('Unexpected index for SPSA.')

		_, deltas = self.get_spsa_params()
		counter = 0
		for param, value in self._method.raw_params.items():
			raw_data_perturbations[param] = value * (1 + direction * deltas[counter])
			counter += 1
		for param, value in self._method.sim_params.items():
			sim_data_perturbations[param] = value * (1 + direction * deltas[counter])
			counter += 1

		return raw_data_perturbations, sim_data_perturbations

	def n_variants_per_iteration(self) -> int:
		"""Number of variants (modified parameter sets) for each iteration"""
		return 2  # high and low change

	### Class specific functions ###

	def get_spsa_params(self):
		it = self.iteration + 1  # prevent 0 for first iteration
		np.random.seed(it)
		at = self.learning_rate / it**self.alpha
		ct = self.parameter_step / it**self.gamma
		deltas = ct * (np.random.rand(self._method.n_parameters) * 2 - 1)

		return at, deltas

