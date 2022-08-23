"""
Adjusts expression of rRNA and rProtein different fixed ppGpp concentrations to
compare which ones are limiting growth at low or high ppGpp concentrations.
Shares a lot of code with ppgpp_limitations

Modifies:
	Attributes from ppgpp_conc variant
	Attributes from adjust_final_expression
	sim_data.ppgpp_ramp: added class to include a ppGpp concentration ramp

Expected variant indices (depends on PPGPP_VARIANTS and FACTORS):
	0: control at low ppGpp concentration
	1-3: adjust rRNA expression at low ppGpp concentration
	4-6: adjust rProtein expression at low ppGpp concentration
	7-9: adjust rRNA and rProtein expression at low ppGpp concentration
	10-12: adjust rRNA (2x rProtein adjustment) and rProtein expression at low ppGpp concentration
	13-15: adjust rRNA (5x rProtein adjustment) and rProtein expression at low ppGpp concentration
	16: control at normal ppGpp concentration
	17-31: adjustments at normal ppGpp concentration
	32: control at high ppGpp concentration
	33-47: adjustments at high ppGpp concentration
"""

import numpy as np

from .ppgpp_conc import ppgpp_conc
from .ppgpp_limitations import ppGpp
from wholecell.utils import units


PPGPP_VARIANTS = [1, 4, 8]  # Corresponds to low, control and high ppGpp from ppgpp_conc variant
FACTORS = [
	[1.1, 1.25, 1.5],
	[1.1, 1.25, 1.5],
	[1.1, 1.25, 1.5],
	[1.1, 1.25, 1.5],
	[1.1, 1.25, 1.5],
]


def adjust_rrna(sim_data, factor):
	rna_data = sim_data.process.transcription.rna_data
	cistron_indices = np.where(rna_data['is_rRNA'])[0]
	factors = [factor] * len(cistron_indices)
	sim_data.adjust_final_expression(cistron_indices, factors)

def adjust_rprotein(sim_data, factor):
	rna_data = sim_data.process.transcription.rna_data
	cistron_indices = np.where(rna_data['includes_ribosomal_protein'])[0]
	factors = [factor] * len(cistron_indices)
	sim_data.adjust_final_expression(cistron_indices, factors)

def adjust_both(sim_data, factor, scale=1):
	rrna_factor = 1 + (factor - 1) * scale
	adjust_rrna(sim_data, rrna_factor)
	adjust_rprotein(sim_data, factor)

def adjust_double(sim_data, factor):
	adjust_both(sim_data, factor, scale=2)

def adjust_5x(sim_data, factor):
	adjust_both(sim_data, factor, scale=5)

ADJUSTMENTS = [
	adjust_rrna,
	adjust_rprotein,
	adjust_both,
	adjust_double,
	adjust_5x,
	]

def split_index(index):
	n_factors = np.unique([len(factors) for factors in FACTORS])
	if len(n_factors) != 1:
		raise ValueError('The number of factors for each adjustment is expected to be the same.')
	n_factors = n_factors[0]
	n_adjustments = len(ADJUSTMENTS)
	n_variants_per_ppgpp = n_factors * n_adjustments + 1  # +1 for control for each ppGpp level

	ppgpp_index = PPGPP_VARIANTS[index // n_variants_per_ppgpp]
	variant_index = index % n_variants_per_ppgpp - 1
	control = variant_index == -1
	adjustment_index = variant_index // n_factors
	factor_index = variant_index % n_factors

	return ppgpp_index, control, adjustment_index, factor_index

def get_adjustment(index):
	_, _, adjustment_index, _ = split_index(index)
	return ADJUSTMENTS[adjustment_index]

def get_factor(index):
	_, control, adjustment_index, factor_index = split_index(index)
	return 1 if control else FACTORS[adjustment_index][factor_index]

def plot_split(index):
	ppgpp_index, control, adjustment_index, factor_index = split_index(index)
	condition = ppgpp_index * len(ADJUSTMENTS) + adjustment_index
	factor = get_factor(index)
	if control:
		factor = 0
	elif factor < 1:
		factor = factor_index - np.sum(np.array(FACTORS[adjustment_index]) < 1)
	else:
		factor = factor_index - np.sum(np.array(FACTORS[adjustment_index]) < 1) + 1
	return condition, factor

def ppgpp_limitations_ribosome(sim_data, index):
	ppgpp_original = sim_data.growth_rate_parameters.get_ppGpp_conc(sim_data.doubling_time)
	ppgpp_index, control, _, _ = split_index(index)

	ppgpp_desc, sim_data = ppgpp_conc(sim_data, ppgpp_index)

	# Get the new concentration after ppgpp_conc updates get_ppGpp_conc function
	ppgpp_rate = 0.00001 * units.mmol / units.L  # per s
	ppgpp_new = sim_data.growth_rate_parameters.get_ppGpp_conc(sim_data.doubling_time)

	# Set a concentration ramp for a gradual change over time for better stability
	sim_data.ppgpp_ramp = ppGpp(ppgpp_original, ppgpp_new, ppgpp_rate)

	if control:
		short_name = 'control'
		desc = 'No parameter changes'
	else:
		adjustment = get_adjustment(index)
		factor = get_factor(index)
		adjustment(sim_data, factor)

		short_name = f'{adjustment.__name__}:{factor}x'
		desc = f'{adjustment.__name__.split("_")} by {factor}x'

	variant_desc = dict(
		shortName=f'{ppgpp_desc["shortName"]} {short_name}',
		desc=f'{desc} with {ppgpp_desc["desc"]}'
		)
	return variant_desc, sim_data
