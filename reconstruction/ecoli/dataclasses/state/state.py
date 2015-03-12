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

from reconstruction.ecoli.dataclasses.state.stateFunctions import createIdsWithCompartments, createMassesByCompartments

import re
import numpy as np

class State(object):
	""" State """

	def __init__(self, raw_data, sim_data):

		self.bulkMolecules = BulkMolecules(raw_data, sim_data)
		self.bulkChromosome = BulkChromosome(raw_data, sim_data)
		self.uniqueMolecules = UniqueMolecules(raw_data, sim_data)

		self._buildBulkMolecules(raw_data, sim_data)
		self._buildBulkChromosome(raw_data, sim_data)
		self._buildUniqueMolecules(raw_data, sim_data)
		self._buildCompartments(raw_data, sim_data)


	def _buildBulkMolecules(self, raw_data, sim_data):

		# Set metabolites
		metaboliteIds = createIdsWithCompartments(raw_data.metabolites)
		metaboliteMasses = units.g / units.mol * createMassesByCompartments(raw_data.metabolites)

		self.bulkMolecules.addToBulkState(metaboliteIds, metaboliteMasses)

		# Set RNA
		rnaIds = createIdsWithCompartments(raw_data.rnas)
		rnaMasses = units.g / units.mol * createMassesByCompartments(raw_data.rnas)

		self.bulkMolecules.addToBulkState(rnaIds, rnaMasses)

		# Set proteins
		proteinIds = createIdsWithCompartments(raw_data.proteins)
		proteinMasses = units.g / units.mol * createMassesByCompartments(raw_data.proteins)

		self.bulkMolecules.addToBulkState(proteinIds, proteinMasses)

		# Set complexes
		complexIds = createIdsWithCompartments(raw_data.proteinComplexes)
		complexMasses = units.g / units.mol * createMassesByCompartments(raw_data.proteinComplexes)

		self.bulkMolecules.addToBulkState(complexIds, complexMasses)

		# Set polymerized
		polymerizedIds = createIdsWithCompartments(raw_data.polymerized)
		polymerizedMasses = units.g / units.mol * createMassesByCompartments(raw_data.polymerized)

		self.bulkMolecules.addToBulkState(polymerizedIds, polymerizedMasses)

	def _buildBulkChromosome(self, raw_data, sim_data):
		# Set genes
		geneIds = [x['id'] for x in raw_data.genes]
		geneMasses = units.g / units.mol * np.zeros((len(geneIds), len(sim_data.molecular_weight_order)), np.float64)

		self.bulkChromosome.addToBulkState(geneIds, geneMasses)

	def _buildUniqueMolecules(self, raw_data, sim_data):
		# Add active RNA polymerase
		rnaPolyComplexMass = self.bulkMolecules.bulkData["mass"][self.bulkMolecules.bulkData["id"] == "APORNAP-CPLX[c]"]
		rnaPolyAttributes = {
				'rnaIndex' : 'i8',
				'transcriptLength' : 'i8'
				}
		self.uniqueMolecules.addToUniqueState('activeRnaPoly', rnaPolyAttributes, rnaPolyComplexMass)

		# Add active ribosome
		# TODO: This is a bad hack that works because in the fitter
		# I have forced expression to be these subunits only
		ribosome30SMass = self.bulkMolecules.bulkData["mass"][
		self.bulkMolecules.bulkData["id"] == sim_data.moleculeGroups.s30_fullComplex[0]
			]
		ribosome50SMass = self.bulkMolecules.bulkData["mass"][
		self.bulkMolecules.bulkData["id"] == sim_data.moleculeGroups.s50_fullComplex[0]
			]
		ribosomeMass = ribosome30SMass + ribosome50SMass
		ribosomeAttributes = {
				'proteinIndex' : 'i8',
				'peptideLength': 'i8'
				}
		self.uniqueMolecules.addToUniqueState('activeRibosome', ribosomeAttributes, ribosomeMass)

		# Add active DNA polymerase
		dnaPolyMass = units.g / units.mol * np.zeros_like(rnaPolyComplexMass) # NOTE: dnaPolymerases currently have no mass
		dnaPolymeraseAttributes = {
				'chromosomeLocation' : 'i8',
				'directionIsPositive' : 'bool',
				'isLeading' : 'bool'
				}
		self.uniqueMolecules.addToUniqueState('dnaPolymerase', dnaPolymeraseAttributes, dnaPolyMass)

	def _buildCompartments(self, raw_data, sim_data):
		compartmentData = np.empty(len(raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in raw_data.compartments]
		self.compartments = compartmentData
