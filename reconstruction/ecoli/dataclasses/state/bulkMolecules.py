"""
SimulationData for bulk molecules state

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class BulkMolecules(object):
	""" BulkMolecules """

	def __init__(self, raw_data, sim_data):
		self._buildBulkMolecules(raw_data, sim_data)


	def _buildBulkMolecules(self, raw_data, sim_data):

		# Instantiate bulkMolecules
		bulkMolecules = np.zeros(
			0,
			dtype = [
				("id", "a50"),
				("mass", "{}f8".format(len(sim_data.molecular_weight_order))),
				]
			)

		# Save compartment abbreviations
		compartmentAbbreviations = np.array([compartment['abbrev'] for compartment in raw_data.compartments])

		# Set metabolites
		metaboliteIds = np.array([metabolite['id'] for metabolite in raw_data.metabolites])
		metaboliteIdsByCompartment = self._createIdsInAllCompartments(metaboliteIds, compartmentAbbreviations)

		metaboliteMassIdxs = np.empty(metaboliteIdsByCompartment.size, np.int64)
		metaboliteMassIdxs.fill(sim_data.molecular_weight_order["metabolite"])
		for index, metabolite in enumerate(metaboliteIdsByCompartment):
			if metabolite.startswith("H2O["):
				metaboliteMassIdxs[index] = sim_data.molecular_weight_order["water"]

		metaboliteMasses = np.zeros((metaboliteIdsByCompartment.size, len(sim_data.molecular_weight_order)), np.float64)
		metaboliteMasses[np.arange(metaboliteIdsByCompartment.size), metaboliteMassIdxs] = [
			metabolite['mw7.2']
			for compartmentIndex in range(compartmentAbbreviations.size)
			for metabolite in raw_data.metabolites
			]

		bulkMolecules = self.addToBulkState(bulkMolecules, metaboliteIdsByCompartment, metaboliteMasses)

		# Set RNA
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location']) for rna in raw_data.rnas]
		rnaMasses = np.array([rna['mw'] for rna in raw_data.rnas])

		bulkMolecules = self.addToBulkState(bulkMolecules, rnaIds, rnaMasses)

		# Set proteins
		# TODO: Change protein masses to be a vector in the flat TSV file like RNA masses
		proteinIds = ['{}[{}]'.format(protein['id'], protein['location']) for protein in raw_data.proteins]
		proteinMasses = np.zeros((len(proteinIds), len(sim_data.molecular_weight_order)), np.float64)
		proteinMasses[np.arange(len(proteinIds)), sim_data.molecular_weight_order["protein"]] = [protein['mw'] for protein in raw_data.proteins]

		bulkMolecules = self.addToBulkState(bulkMolecules, proteinIds, proteinMasses)

		# Set complexes
		complexIds = ['{}[{}]'.format(complex_['id'],complex_['location']) for complex_ in raw_data.proteinComplexes]
		complexMasses = np.array([complex_['mw'] for complex_ in raw_data.proteinComplexes])

		bulkMolecules = self.addToBulkState(bulkMolecules, complexIds, complexMasses)
		
		# Set polymerized
		polymerizedIDs = [entry["id"] for entry in raw_data.polymerized]

		polymerizedIDsByCompartment = [
			'{}[{}]'.format(polymerizedID, compartmentAbbreviation)
			for compartmentAbbreviation in compartmentAbbreviations
			for polymerizedID in polymerizedIDs
			]

		polymerizedMasses = np.zeros((len(polymerizedIDsByCompartment), len(sim_data.molecular_weight_order)), np.int64)
		polymerizedMassesIndexes = [
			entry["mass key"]
			for compartmentAbbreviation in compartmentAbbreviations
			for entry in raw_data.polymerized
			]
		polymerizedMasses[np.arange(len(polymerizedIDsByCompartment)), polymerizedMassesIndexes] = [
			entry["mw"]
			for compartmentAbbreviation in compartmentAbbreviations
			for entry in raw_data.polymerized
			]

		bulkMolecules = self.addToBulkState(bulkMolecules, polymerizedIDsByCompartment, polymerizedMasses)

		# Add units to values
		field_units = {
			"id"		:	None,
			"mass"				:	units.g / units.mol,
			}

		self.bulkMolecules = UnitStructArray(bulkMolecules, field_units)


	## Helper Functions ##
	def addToBulkState(self, bulkState, ids, masses):
		newAddition = np.zeros(
			len(ids),
			dtype = [
				("id", "a50"),
				("mass", "{}f8".format(masses.shape[1])), # TODO: Make this better
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
