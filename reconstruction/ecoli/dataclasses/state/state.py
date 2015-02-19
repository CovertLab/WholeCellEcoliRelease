"""
SimulationData state associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.bulkMolecules import BulkMolecules
from reconstruction.ecoli.dataclasses.state.bulkChromosome import BulkChromosome
from reconstruction.ecoli.dataclasses.state.uniqueMolecules import UniqueMolecules

from reconstruction.ecoli.dataclasses.state.bulkStateFunctions import addToBulkState, createIdsInAllCompartments

import re
import numpy as np

class State(object):
	""" State """

	def __init__(self, raw_data, sim_data):
		self.addToBulkState = addToBulkState
		self.createIdsInAllCompartments = createIdsInAllCompartments

		self.bulkMolecules = BulkMolecules(raw_data, sim_data)
		self.bulkChromosome = BulkChromosome(raw_data, sim_data)
		#self.uniqueMolecules = UniqueMolecules(raw_data, sim_data)

		self._buildBulkMolecules(raw_data, sim_data)
		self._buildBulkChromosome(raw_data, sim_data)
		#self._buildUniqueMolecules(raw_data, sim_data)
		self._buildCompartments(raw_data, sim_data)


	def _buildBulkMolecules(self, raw_data, sim_data):
		# Save compartment abbreviations
		compartmentAbbreviations = np.array([compartment['abbrev'] for compartment in raw_data.compartments])

		# Set metabolites
		metaboliteIds = np.array([metabolite['id'] for metabolite in raw_data.metabolites])
		metaboliteIdsByCompartment = createIdsInAllCompartments(metaboliteIds, compartmentAbbreviations)

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

		self.bulkMolecules.bulkMoleculesData = addToBulkState(self.bulkMolecules.bulkMoleculesData, metaboliteIdsByCompartment, metaboliteMasses)

		# Set RNA
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location']) for rna in raw_data.rnas]
		rnaMasses = np.array([rna['mw'] for rna in raw_data.rnas])

		self.bulkMolecules.bulkMoleculesData = addToBulkState(self.bulkMolecules.bulkMoleculesData, rnaIds, rnaMasses)

		# Set proteins
		# TODO: Change protein masses to be a vector in the flat TSV file like RNA masses
		proteinIds = ['{}[{}]'.format(protein['id'], protein['location']) for protein in raw_data.proteins]
		proteinMasses = np.zeros((len(proteinIds), len(sim_data.molecular_weight_order)), np.float64)
		proteinMasses[np.arange(len(proteinIds)), sim_data.molecular_weight_order["protein"]] = [protein['mw'] for protein in raw_data.proteins]

		self.bulkMolecules.bulkMoleculesData = addToBulkState(self.bulkMolecules.bulkMoleculesData, proteinIds, proteinMasses)

		# Set complexes
		complexIds = ['{}[{}]'.format(complex_['id'],complex_['location']) for complex_ in raw_data.proteinComplexes]
		complexMasses = np.array([complex_['mw'] for complex_ in raw_data.proteinComplexes])

		self.bulkMolecules.bulkMoleculesData = addToBulkState(self.bulkMolecules.bulkMoleculesData, complexIds, complexMasses)
		
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

		self.bulkMolecules.bulkMoleculesData = addToBulkState(self.bulkMolecules.bulkMoleculesData, polymerizedIDsByCompartment, polymerizedMasses)

		# Add units to values
		field_units = {
			"id"		:	None,
			"mass"				:	units.g / units.mol,
			}

		self.bulkMolecules.bulkMoleculesData = UnitStructArray(self.bulkMolecules.bulkMoleculesData, field_units)

	def _buildBulkChromosome(self, raw_data, sim_data):
		pass

	def _buildUniqueMolecules(self, raw_data, sim_data):
		pass

	def _buildCompartments(self, raw_data, sim_data):
		compartmentData = np.empty(len(raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in raw_data.compartments]
		self.compartments = compartmentData
