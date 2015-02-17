"""
SimulationData state associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

import reconstruction.ecoli.dataclasses.dataclass
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

import re
import numpy as np

class State(reconstruction.ecoli.dataclasses.dataclass.DataClass):
	""" State """

	def __init__(self, simData):
		super(State, self).__init__(simData)

		self._buildCompartments()
		self._buildBulkMolecules()
		self._buildBulkChromosome()


	def _buildCompartments(self):
		compartmentData = np.empty(len(self._simData.raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in self._simData.raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in self._simData.raw_data.compartments]
		self.compartments = compartmentData
		self.n_compartments = compartmentData.size


	def _addToBulkState(self, bulkState, ids, masses):
		newAddition = np.zeros(
			len(ids),
			dtype = [
				("id", "a50"),
				("mass", "{}f8".format(len(self._simData.molecular_weight_order))),
				]
			)

		newAddition["id"] = ids
		newAddition["mass"] = masses
		return np.hstack((bulkState, newAddition))

	def _createIdsInAllCompartments(self, ids, compartments):
		idsByCompartment = [
			'{}[{}]'.format(i, c)
			for c in compartments
			for i in ids
			]
		return np.array(idsByCompartment)

	def _buildBulkMolecules(self):

		# Instantiate bulkMolecules
		bulkMolecules = np.zeros(
			0,
			dtype = [
				("id", "a50"),
				("mass", "{}f8".format(len(self._simData.molecular_weight_order))),
				]
			)

		# Save compartment abbreviations
		compartmentAbbreviations = np.array([compartment['abbrev'] for compartment in self._simData.raw_data.compartments])

		# Set metabolites
		metaboliteIds = np.array([metabolite['id'] for metabolite in self._simData.raw_data.metabolites])
		metaboliteIdsByCompartment = self._createIdsInAllCompartments(metaboliteIds, compartmentAbbreviations)

		metaboliteMassIdxs = np.empty(metaboliteIdsByCompartment.size, np.int64)
		metaboliteMassIdxs.fill(self._simData.molecular_weight_order["metabolite"])
		for index, metabolite in enumerate(metaboliteIdsByCompartment):
			if metabolite.startswith("H2O["):
				metaboliteMassIdxs[index] = self._simData.molecular_weight_order["water"]

		metaboliteMasses = np.zeros((metaboliteIdsByCompartment.size, len(self._simData.molecular_weight_order)), np.float64)
		metaboliteMasses[np.arange(metaboliteIdsByCompartment.size), metaboliteMassIdxs] = [
			metabolite['mw7.2']
			for compartmentIndex in range(compartmentAbbreviations.size)
			for metabolite in self._simData.raw_data.metabolites
			]

		bulkMolecules = self._addToBulkState(bulkMolecules, metaboliteIdsByCompartment, metaboliteMasses)

		# Set RNA
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location']) for rna in self._simData.raw_data.rnas]
		rnaMasses = [rna['mw'] for rna in self._simData.raw_data.rnas]

		bulkMolecules = self._addToBulkState(bulkMolecules, rnaIds, rnaMasses)

		# Set proteins
		# TODO: Change protein masses to be a vector in the flat TSV file like RNA masses
		proteinIds = ['{}[{}]'.format(protein['id'], protein['location']) for protein in self._simData.raw_data.proteins]
		proteinMasses = np.zeros((len(proteinIds), len(self._simData.molecular_weight_order)), np.float64)
		proteinMasses[np.arange(len(proteinIds)), self._simData.molecular_weight_order["protein"]] = [protein['mw'] for protein in self._simData.raw_data.proteins]

		bulkMolecules = self._addToBulkState(bulkMolecules, proteinIds, proteinMasses)

		# Set complexes
		complexIds = ['{}[{}]'.format(complex_['id'],complex_['location']) for complex_ in self._simData.raw_data.proteinComplexes]
		complexMasses = [complex_['mw'] for complex_ in self._simData.raw_data.proteinComplexes]

		bulkMolecules = self._addToBulkState(bulkMolecules, complexIds, complexMasses)
		
		# Set polymerized
		polymerizedIDs = [entry["id"] for entry in self._simData.raw_data.polymerized]

		polymerizedIDsByCompartment = [
			'{}[{}]'.format(polymerizedID, compartmentAbbreviation)
			for compartmentAbbreviation in compartmentAbbreviations
			for polymerizedID in polymerizedIDs
			]

		polymerizedMasses = np.zeros((len(polymerizedIDsByCompartment), len(self._simData.molecular_weight_order)), np.int64)
		polymerizedMassesIndexes = [
			entry["mass key"]
			for compartmentAbbreviation in compartmentAbbreviations
			for entry in self._simData.raw_data.polymerized
			]
		polymerizedMasses[np.arange(len(polymerizedIDsByCompartment)), polymerizedMassesIndexes] = [
			entry["mw"]
			for compartmentAbbreviation in compartmentAbbreviations
			for entry in self._simData.raw_data.polymerized
			]

		bulkMolecules = self._addToBulkState(bulkMolecules, polymerizedIDsByCompartment, polymerizedMasses)

		# Add units to values
		field_units = {
			"id"		:	None,
			"mass"				:	units.g / units.mol,
			}

		self.bulkMolecules = UnitStructArray(bulkMolecules, field_units)


	def _buildBulkChromosome(self):
		bulkChromosome = np.zeros(0,
			dtype = [("id", 			"a50"),
					("mass", "{}f8".format(len(self._simData.molecular_weight_order)))
					]
					)

		# Set genes
		geneIds = [x['id'] for x in self._simData.raw_data.genes]
		geneMasses = np.zeros((len(geneIds), len(self._simData.molecular_weight_order)), np.float64)

		bulkChromosome = self._addToBulkState(bulkChromosome, geneIds, geneMasses)

		# Add units to values
		field_units = {
			"id"			:	None,
			"mass"					:	units.g / units.mol,
			}
		self.bulkChromosome = UnitStructArray(bulkChromosome, field_units)
