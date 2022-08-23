"""
Simulation data manual adjustments necessary for proper metabolism.

Raw tsv files in 'reconstruction/ecoli/flat/adjustments'

"""

from __future__ import absolute_import, division, print_function

class Adjustments(object):

	def __init__(self, raw_data):
		self._add_adjustments(raw_data)

	def _add_adjustments(self, raw_data):
		self.translation_efficiencies_adjustments = {
			adj["name"]: adj["value"]
			for adj in raw_data.adjustments.translation_efficiencies_adjustments
		}
		self.rna_expression_adjustments = {
			adj["name"]: adj["value"]
			for adj in raw_data.adjustments.rna_expression_adjustments
		}
		self.rna_deg_rates_adjustments = {
			adj["name"]: adj["value"]
			for adj in raw_data.adjustments.rna_deg_rates_adjustments
		}
		self.protein_deg_rates_adjustments = {
			adj["name"]: (adj["value"] if not isinstance(adj["value"], str) else eval(adj["value"]))  # eval fractions
			for adj in raw_data.adjustments.protein_deg_rates_adjustments
		}
		self.relative_metabolite_concentrations_changes = {
			adj["media"]: {adj["metabolite"]: adj["fold_change"]}
			for adj in raw_data.adjustments.relative_metabolite_concentrations_changes
		}
