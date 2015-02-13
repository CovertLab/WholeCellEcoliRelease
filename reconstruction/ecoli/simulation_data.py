"""
SimulationData for Ecoli

Raw data processed into forms convienent for whole-cell modeling

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""
from __future__ import division

import numpy as np

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

# Data classes
from reconstruction.ecoli.dataclasses.getterFunctions import getterFunctions
from reconstruction.ecoli.dataclasses.moleculeGroups import moleculeGroups
from reconstruction.ecoli.dataclasses.replication import Replication
from reconstruction.ecoli.dataclasses.state import State

# Unit imports
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class SimulationDataEcoli(object):
	""" SimulationDataEcoli """

	def __init__(self):
		# Raw data
		self.raw_data = KnowledgeBaseEcoli()

		# Data classes
		self.getter = getterFunctions(self)
		self.moleculeGroups = moleculeGroups(self)
		self.replication = Replication(self)
		self.state = State(self)

		# Functions
		self._calculateUsefulParameters()
		self._buildAllMasses()
		
	def _calculateUsefulParameters(self):
		self.sizeBulkMolecules = len(self.raw_data.rnas) + len(self.raw_data.proteins) + len(self.raw_data.proteinComplexes) + len(self.raw_data.metabolites) + len(self.raw_data.polymerized)

	def _buildAllMasses(self):
		allMass = np.empty(self.sizeBulkMolecules,
			dtype = [
					('id',		'a50'),
					('mass',	"f8")
					]
			)

		listMass = []
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in self.raw_data.rnas])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in self.raw_data.proteins])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in self.raw_data.proteinComplexes])
		listMass.extend([(x['id'],np.sum(x['mw7.2'])) for x in self.raw_data.metabolites])
		listMass.extend([(x['id'],np.sum(x['mw'])) for x in self.raw_data.polymerized])

		allMass[:] = listMass

		field_units = {
			'id'		:	None,
			'mass'		:	units.g / units.mol,
			}

		self._allMass = UnitStructArray(allMass, field_units)
