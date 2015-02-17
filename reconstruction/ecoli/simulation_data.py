"""
SimulationData for Ecoli

Raw data processed into forms convienent for whole-cell modeling

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""
from __future__ import division

import numpy as np
import collections

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

# Data classes
from reconstruction.ecoli.dataclasses.getterFunctions import getterFunctions
from reconstruction.ecoli.dataclasses.moleculeGroups import moleculeGroups
from reconstruction.ecoli.dataclasses.state.state import State
from reconstruction.ecoli.dataclasses.process.process import Process

# Unit imports
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class SimulationDataEcoli(object):
	""" SimulationDataEcoli """

	def __init__(self):
		# Raw data
		raw_data = KnowledgeBaseEcoli()
		self._addHardCodedAttributes()

		# Data classes
		self.getter = getterFunctions(raw_data, self)
		self.moleculeGroups = moleculeGroups(raw_data, self)
		self.process = Process(raw_data, self)
		self.state = State(raw_data, self)

		# Functions
		self._buildAllMasses(raw_data)

	def _buildAllMasses(self, raw_data):
		size = len(raw_data.rnas) + len(raw_data.proteins) + len(raw_data.proteinComplexes) + len(raw_data.metabolites) + len(raw_data.polymerized)
		allMass = np.empty(size,
			dtype = [
					('id',		'a50'),
					('mass',	"f8")
					]
			)

		listMass = []
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in raw_data.rnas])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in raw_data.proteins])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in raw_data.proteinComplexes])
		listMass.extend([(x['id'],np.sum(x['mw7.2'])) for x in raw_data.metabolites])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in raw_data.polymerized])

		allMass[:] = listMass

		field_units = {
			'id'		:	None,
			'mass'		:	units.g / units.mol,
			}

		self._allMass = UnitStructArray(allMass, field_units)

	def _addHardCodedAttributes(self):
		self.molecular_weight_keys = [
			'23srRNA',
			'16srRNA',
			'5srRNA',
			'tRNA',
			'mRNA',
			'miscRNA',
			'protein',
			'metabolite',
			'water',
			'DNA',
			'RNA' # nonspecific RNA
			]

		self.molecular_weight_order = collections.OrderedDict([
			(key, index) for index, key in enumerate(self.molecular_weight_keys)
			])