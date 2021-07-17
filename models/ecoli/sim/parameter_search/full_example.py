"""
An example of using the parameter search framework to optimize over two
parameters, one modified before the parca and one modified after.  Each
iteration will need to run the parca because _raw_params are specified.
"""

import os
import pickle

from models.ecoli.sim.parameter_search.base_parameter_search import BaseParameterSearch, RawParameter, SimParameter


class FullExample(BaseParameterSearch):
	_raw_params = (RawParameter('metabolite_concentrations', {'Metabolite': 'TRP'}, ['Park Concentration', 'Lempp Concentration'], 'trp conc'),)
	_sim_params = (SimParameter('constants.test'),)
	_init_sim_params = {'constants.test': 1}
	sims_to_run = (
		{
			'jit': False,
			'length_sec': 2,
			'variant': ('condition', 0),
		},
	)

	def get_objective(self, sim_out_dirs, sim_data_files):
		objectives = []
		for sim_dir, sim_data_file in zip(sim_out_dirs, sim_data_files):
			out_dir = os.path.join(sim_dir, '000000', 'generation_000000', '000000', 'simOut')  # TODO: function to get all simOut dirs
			with open(sim_data_file, 'rb') as f:
				sim_data = pickle.load(f)

			aa_ids = sim_data.molecule_groups.amino_acids
			trp_idx = aa_ids.index('TRP[c]')
			trp_conc = self.read_column(out_dir, 'GrowthLimits', 'aa_supply_aa_conc')[1:, trp_idx].mean()

			objective = sim_data.constants.test**2 + (trp_conc - 0.05)**2
			objectives.append(objective)

		return objectives
