"""
SimulationData for metabolism process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
import numpy as np

class Metabolism(object):
	""" Metabolism """

	def __init__(self, raw_data, sim_data):
		self._buildBiomass(raw_data, sim_data)

	def _buildBiomass(self, raw_data, sim_data):
		coreBiomassIds = []
		coreBiomassConcentration = []
		for metabolite in raw_data.metabolites:
			for location in metabolite['core_location']:
				coreBiomassIds.append('{}[{}]'.format(metabolite['id'], location))
			for concentration in metabolite['core_concentration']:
				coreBiomassConcentration.append(concentration)

		wildtypeBiomassIds = []
		wildtypeBiomassConcentration = []
		for metabolite in raw_data.metabolites:
			for location in metabolite['wildtype_location']:
				wildtypeBiomassIds.append('{}[{}]'.format(metabolite['id'], location))
			for concentration in metabolite['wildtype_concentration']:
				wildtypeBiomassConcentration.append(concentration)
				
		coreBiomassData = np.zeros(len(coreBiomassIds), dtype = [('id', 'a50'), ('biomassFlux', 'f8')])
		wildtypeBiomassData = np.zeros(len(wildtypeBiomassIds), dtype = [('id', 'a50'), ('biomassFlux',		'f8')])

		coreBiomassData['id'] = coreBiomassIds
		coreBiomassData['biomassFlux'] = coreBiomassConcentration
		wildtypeBiomassData['id'] = wildtypeBiomassIds
		wildtypeBiomassData['biomassFlux'] = wildtypeBiomassConcentration

		field_units = {'id' : None,
				'biomassFlux' : units.mmol / units.g}
		self.coreBiomass 		= UnitStructArray(coreBiomassData, field_units)
		self.wildtypeBiomass 	= UnitStructArray(wildtypeBiomassData, field_units)

