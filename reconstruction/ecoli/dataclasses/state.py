"""
SimulationData state associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

import reconstruction.ecoli.dataclasses.dataclass

import re
import numpy as np

class State(reconstruction.ecoli.dataclasses.dataclass.DataClass):
	""" State """

	def __init__(self, simData):
		super(State, self).__init__(simData)

		self._buildCompartments()


	def _buildCompartments(self):
		compartmentData = np.empty(len(self._simData.raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in self._simData.raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in self._simData.raw_data.compartments]
		self.compartments = compartmentData
		self.n_compartments = compartmentData.size

	def _buildBulkMolecules(self):
		# TODO: modularize this logic

		size = (
			len(self._metabolites)*len(self._compartmentList)
			+ len(self._rnas)
			+ len(self._proteins)
			+ len(self._proteinComplexes)
			+ len(self._polymerized)*len(self._compartmentList)
			)

		bulkMolecules = np.zeros(
			size,
			dtype = [
				("moleculeId", "a50"),
				('compartment',	 "a1"),
				("mass", "{}f8".format(len(MOLECULAR_WEIGHT_ORDER))),
				]
			)

		# Set metabolites
		lastMetaboliteIdx = len(self._metabolites) * len(self._compartmentList)

		compartmentAbbreviations = [compartment['abbrev'] for compartment in self._compartmentList]
		metaboliteIds = [metabolite['id'] for metabolite in self._metabolites]

		bulkMolecules['moleculeId'][0:lastMetaboliteIdx] = [
			'{}[{}]'.format(metaboliteId, compartmentAbbreviation)
			for compartmentAbbreviation in compartmentAbbreviations
			for metaboliteId in metaboliteIds
			]

		metaboliteMassIdxs = np.empty(lastMetaboliteIdx, np.int64)

		metaboliteMassIdxs.fill(MOLECULAR_WEIGHT_ORDER["metabolite"])

		for index, metabolite in enumerate(bulkMolecules[0:lastMetaboliteIdx]):
			if metabolite["moleculeId"].startswith("H2O["):
				metaboliteMassIdxs[index] = MOLECULAR_WEIGHT_ORDER["water"]

		bulkMolecules['mass'][np.arange(lastMetaboliteIdx), metaboliteMassIdxs] = [
			metabolite['mw7.2']
			for compartmentIndex in range(len(self._compartmentList))
			for metabolite in self._metabolites
			]

		# Set RNA
		lastRnaIdx = len(self._rnas) + lastMetaboliteIdx

		bulkMolecules['moleculeId'][lastMetaboliteIdx:lastRnaIdx] = [
			'{}[{}]'.format(rna['id'], rna['location']) for rna in self._rnas
			]

		bulkMolecules['mass'][lastMetaboliteIdx:lastRnaIdx, :] = [
			rna['mw'] for rna in self._rnas
			]

		# Set proteins
		lastProteinMonomerIdx = len(self._proteins) + lastRnaIdx

		bulkMolecules['moleculeId'][lastRnaIdx:lastProteinMonomerIdx] = [
			'{}[{}]'.format(protein['id'], protein['location'])
			for protein in self._proteins
			]

		bulkMolecules['mass'][lastRnaIdx:lastProteinMonomerIdx, MOLECULAR_WEIGHT_ORDER["protein"]] = [
			protein['mw'] for protein in self._proteins
			]

		# Set complexes
		lastComplexIdx = len(self._proteinComplexes) + lastProteinMonomerIdx

		bulkMolecules['moleculeId'][lastProteinMonomerIdx:lastComplexIdx] = [
			'{}[{}]'.format(complex_['id'],complex_['location']) for complex_ in self._proteinComplexes
			]

		bulkMolecules['mass'][lastProteinMonomerIdx:lastComplexIdx, :] = [
			complex_['mw'] for complex_ in self._proteinComplexes
			]

		# Set polymerized

		lastPolymerizedIndex = len(self._polymerized)*len(self._compartmentList) + lastComplexIdx

		polymerizedIDs = [entry["id"] for entry in self._polymerized]

		bulkMolecules["moleculeId"][lastComplexIdx:lastPolymerizedIndex] = [
			'{}[{}]'.format(polymerizedID, compartmentAbbreviation)
			for compartmentAbbreviation in compartmentAbbreviations
			for polymerizedID in polymerizedIDs
			]

		masses = [
			entry["mw"]
			for compartmentAbbreviation in compartmentAbbreviations
			for entry in self._polymerized
			]

		massIndexes = [
			entry["mass key"]
			for compartmentAbbreviation in compartmentAbbreviations
			for entry in self._polymerized
			]

		bulkMolecules["mass"][range(lastComplexIdx, lastPolymerizedIndex), massIndexes] = masses
		# NOTE: the use of range above is intentional

		# Add units to values
		field_units = {
			"moleculeId"		:	None,
			"mass"				:	units.g / units.mol,
			'compartment'		:	None,
			}

		self.bulkMolecules = UnitStructArray(bulkMolecules, field_units)