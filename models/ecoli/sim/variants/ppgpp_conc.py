"""
Set the concentration of ppGpp to be constant during simulations at different
concentrations.

Modifies:
	attributes modified by the condition variant
	sim_data.growth_rate_parameters.get_ppGpp_conc
	sim_data.constants.k_RelA_ppGpp_synthesis
	sim_data.constants.k_SpoT_ppGpp_degradation
	sim_data.constants.k_SpoT_ppGpp_synthesis
	sim_data.process.metabolism.force_constant_ppgpp

Expected variant indices (dependent on FACTORS and CONDITIONS):
	0-3: lower concentrations of ppGpp for first condition
	4: control for first condtiion
	5-9: higher concentrations of ppGpp for first condition
	10-19: ppGpp concentrations for second condition
"""

from .condition import condition
from wholecell.utils import units


FACTORS = [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2]
CONDITIONS = [0, 2]  # minimal and minimal+AA
BASE_FACTOR = [1, 0.4]  # factors that correspond to expected concentration in CONDITIONS


class ppGpp():
	def __init__(self, conc):
		self.conc = conc

	def get_ppGpp_conc(self, doubling_time):
		"""
		Function to replace get_ppGpp_conc in growth_rate_parameters.  Needs to
		have the same function signature but we want to fix the output
		conentration to a new value.
		"""
		return self.conc


def split_index(index):
	n_factors = len(FACTORS)
	condition_index = index // n_factors
	factor_index = index % n_factors
	return CONDITIONS[condition_index], FACTORS[factor_index]

def ppgpp_conc(sim_data, index):
	condition_index, factor = split_index(index)

	_, sim_data = condition(sim_data, condition_index)
	control_conc = sim_data.growth_rate_parameters.get_ppGpp_conc(sim_data.doubling_time)
	new_conc = factor * control_conc

	# Replace function to only return the desired concentration
	sim_data.growth_rate_parameters.get_ppGpp_conc = ppGpp(new_conc).get_ppGpp_conc

	# Set rate constants to 0 to prevent ppGpp concentration changes during sim
	sim_data.constants.k_RelA_ppGpp_synthesis *= 0
	sim_data.constants.k_SpoT_ppGpp_degradation *= 0
	sim_data.constants.k_SpoT_ppGpp_synthesis *= 0

	# Force constant ppGpp concentration in metabolism (added to homeostatic objective)
	sim_data.process.metabolism.force_constant_ppgpp = True

	return dict(
		shortName=f'{sim_data.condition}:ppGpp:{factor}x',
		desc=f'ppGpp conc adjusted by {factor}x to {new_conc.asNumber(units.umol / units.L)} uM in {sim_data.condition}.'
		), sim_data
