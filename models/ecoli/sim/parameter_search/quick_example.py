"""
An example of using the parameter search framework to optimize over one
parameter that is modified after running the parca.  The parca will not need
to be run with each iteration and no sims are specified to run so evaulations
can occur more quickly to test algorithms.
"""

import pickle

from models.ecoli.sim.parameter_search.base_parameter_search import BaseParameterSearch, SimParameter


class QuickExample(BaseParameterSearch):
	_sim_params = (SimParameter('constants.test'),)
	_init_sim_params = {'constants.test': 1}

	def get_objective(self, sim_out_dirs, sim_data_files):
		objectives = []
		for sim_data_file in sim_data_files:
			with open(sim_data_file, 'rb') as f:
				sim_data = pickle.load(f)

			objective = sim_data.constants.test**2
			objectives.append(objective)

		return objectives
