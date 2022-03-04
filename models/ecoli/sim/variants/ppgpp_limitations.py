"""
Adjusts amino acid synthesis parameters and enzymes as well as ribosomes at
different fixed ppGpp concentrations to compare which ones are limiting growth
at low or high ppGpp concentrations.

Modifies:
	Attributes from ppgpp_conc variant
	Attributes from adjust_final_expression
	sim_data.ppgpp_ramp: added class to include a ppGpp concentration ramp
	sim_data.process.metabolism.aa_kcats_fwd
	sim_data.process.metabolism.aa_kis

Expected variant indices (depends on PPGPP_VARIANTS and FACTORS):
	0: control at low ppGpp concentration
	1-9: adjust amino acid KIs at low ppGpp concentration
	10-18: adjust amino acid kcats at low ppGpp concentration
	19-27: adjust amino acid enzyme expression at low ppGpp concentration
	28-36: adjust ribosome expression at low ppGpp concentration
	37: control at normal ppGpp concentration
	38-73: adjustments at normal ppGpp concentration
	74: control at high ppGpp concentration
	75-110: adjustments at high ppGpp concentration
"""

import numpy as np

from .ppgpp_conc import ppgpp_conc
from wholecell.utils import units


PPGPP_VARIANTS = [1, 4, 8]  # Corresponds to low, control and high ppGpp from ppgpp_conc variant
FACTORS = [
	[0.001, 0.01, 0.1, 0.5, 2, 10, 100, 1000, np.inf],
	[0.96, 0.97, 0.98, 0.99, 1.01, 1.02, 1.03, 1.04, 1.05],
	[0.5, 0.75, 0.9, 0.95, 1.05, 1.1, 1.25, 1.5, 2],
	[0.5, 0.75, 0.9, 0.95, 1.05, 1.1, 1.25, 1.5, 2],
]


class ppGpp():
	def __init__(self, original_conc, new_conc, rate):
		self.original_conc = original_conc
		self.new_conc = new_conc
		self.rate = rate

	def time(self):
		return 0

	def get_ppGpp_conc(self, doubling_time):
		"""
		Function to replace get_ppGpp_conc in growth_rate_parameters.  Needs to
		have the same function signature but we want to fix the output
		conentration to a new value.
		"""
		# TODO (travis): handle this error from sim side so that self.time()
		# never raises an error even if called during process initialization
		try:
			t = self.time()
		except Exception:
			t = 0

		if self.new_conc > self.original_conc:
			conc = min(self.new_conc, self.original_conc + self.rate * t)
		else:
			conc = max(self.new_conc, self.original_conc - self.rate * t)
		return conc

	def set_time(self, time_fun):
		self.time = time_fun

def adjust_amino_acids(sim_data, factor):
	sim_data.process.metabolism.aa_kcats_fwd *= factor

def adjust_amino_acid_kis(sim_data, factor):
	sim_data.process.metabolism.aa_kis *= factor

def adjust_enzymes(sim_data, factor):
	metabolism = sim_data.process.metabolism
	complexation = sim_data.process.complexation
	translation = sim_data.process.translation
	transcription = sim_data.process.transcription

	# Enzymes involved in mechanistic amino acid synthesis
	synthesis_enzymes = metabolism.aa_enzymes[metabolism.enzyme_to_amino_acid_fwd.sum(1).astype(bool)]
	synthesis_monomers = sorted({
		subunit
		for enzyme in synthesis_enzymes
		for subunit in complexation.get_monomers(enzyme)['subunitIds']
		})

	# Map monomers to RNA for a knockout
	cistron_to_index = {cistron['id']: i for i, cistron in enumerate(transcription.cistron_data)}
	monomer_to_index = {monomer['id']: cistron_to_index[monomer['cistron_id']] for monomer in translation.monomer_data}
	cistron_indices = [
		monomer_to_index[monomer]
		for monomer in synthesis_monomers
		if monomer in monomer_to_index
		]

	factors = [factor] * len(cistron_indices)
	sim_data.adjust_final_expression(cistron_indices, factors)

def adjust_ribosomes(sim_data, factor):
	rna_data = sim_data.process.transcription.rna_data
	cistron_indices = np.where(rna_data['is_rRNA'] | rna_data['includes_ribosomal_protein'])[0]
	factors = [factor] * len(cistron_indices)
	sim_data.adjust_final_expression(cistron_indices, factors)

ADJUSTMENTS = [
	adjust_amino_acid_kis,
	adjust_amino_acids,
	adjust_enzymes,
	adjust_ribosomes,
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

def ppgpp_limitations(sim_data, index):
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
