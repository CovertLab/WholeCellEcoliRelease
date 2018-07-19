"""
SimulationData state associated data

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.bulkMolecules import BulkMolecules
from reconstruction.ecoli.dataclasses.state.uniqueMolecules import UniqueMolecules

from reconstruction.ecoli.dataclasses.state import stateFunctions as sf

import re
import numpy as np

class InternalState(object):
	""" Internal State """

	def __init__(self, raw_data, sim_data):

		self.bulkMolecules = BulkMolecules(raw_data, sim_data)
		self.uniqueMolecules = UniqueMolecules(raw_data, sim_data)

		self._buildBulkMolecules(raw_data, sim_data)
		self._buildUniqueMolecules(raw_data, sim_data)
		self._buildCompartments(raw_data, sim_data)


	def _buildBulkMolecules(self, raw_data, sim_data):

		# Set metabolites
		metaboliteIds = sf.createIdsWithCompartments(raw_data.metabolites)
		metaboliteMasses = units.g / units.mol * sf.createMetaboliteMassesByCompartments(raw_data.metabolites, 7, 11)

		self.bulkMolecules.addToBulkState(metaboliteIds, metaboliteMasses)

		# Set water
		waterIds = sf.createIdsWithCompartments(raw_data.water)
		waterMasses = units.g / units.mol * sf.createMetaboliteMassesByCompartments(raw_data.water, 8, 11)

		self.bulkMolecules.addToBulkState(waterIds, waterMasses)

		# Set RNA
		rnaIds = sf.createIdsWithCompartments(raw_data.rnas)
		rnaMasses = units.g / units.mol * sf.createMassesByCompartments(raw_data.rnas)

		self.bulkMolecules.addToBulkState(rnaIds, rnaMasses)

		# Set proteins
		proteinIds = sf.createIdsWithCompartments(raw_data.proteins)
		proteinMasses = units.g / units.mol * sf.createMassesByCompartments(raw_data.proteins)

		self.bulkMolecules.addToBulkState(proteinIds, proteinMasses)

		# Set complexes
		complexIds = sf.createIdsWithCompartments(raw_data.proteinComplexes)
		complexMasses = units.g / units.mol * sf.createMassesByCompartments(raw_data.proteinComplexes)

		self.bulkMolecules.addToBulkState(complexIds, complexMasses)

		# Set modified forms
		modifiedFormIds = sf.createIdsWithCompartments(raw_data.modifiedForms)
		modifiedFormMasses = units.g / units.mol * sf.createModifiedFormMassesByCompartments(raw_data.modifiedForms)

		self.bulkMolecules.addToBulkState(modifiedFormIds, modifiedFormMasses)

		# Set chromosome
		chromosomeIds = sf.createIdsWithCompartments(raw_data.chromosome)
		chromosomeMasses = units.g / units.mol * sf.createMassesByCompartments(raw_data.chromosome)

		self.bulkMolecules.addToBulkState(chromosomeIds, chromosomeMasses)

		# Set fragments
		test = []
		for x in raw_data.polymerized:
			if x['is_ntp']:
				if not x['is_end']:
					temp = x
					temp['id'] = x['id'].replace('Polymerized','Fragment')
					test.append(temp)
		fragmentsIds = sf.createIdsWithCompartments(test)
		fragmentsMasses = units.g / units.mol * sf.createMassesByCompartments(test)

		self.bulkMolecules.addToBulkState(fragmentsIds, fragmentsMasses)

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
		self.bulkMolecules.bulkData["id"] == sim_data.moleculeIds.s30_fullComplex
			]
		ribosome50SMass = self.bulkMolecules.bulkData["mass"][
		self.bulkMolecules.bulkData["id"] == sim_data.moleculeIds.s50_fullComplex
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
				'sequenceIdx' : 'i8',
				'sequenceLength' : 'i8',
				'replicationRound' : 'i8',
				'chromosomeIndex' : 'i8',
				}
		self.uniqueMolecules.addToUniqueState('dnaPolymerase', dnaPolymeraseAttributes, dnaPolyMass)

		# Origin of replication
		originMass = units.g / units.mol * np.zeros_like(rnaPolyComplexMass) # NOTE: origins currently have no mass
		originAttributes = {
				'chromosomeIndex': 'i8',
				}
		self.uniqueMolecules.addToUniqueState('originOfReplication', originAttributes, originMass)

		# Full chromosome
		fullChromosomeMass = units.g / units.mol * np.zeros_like(rnaPolyComplexMass) # NOTE: origins currently have no mass
		fullChromosomeAttributes = {"division_time" : "f8"}
		self.uniqueMolecules.addToUniqueState('fullChromosome', fullChromosomeAttributes, fullChromosomeMass)


	def _buildCompartments(self, raw_data, sim_data):
		compartmentData = np.empty(len(raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in raw_data.compartments]
		self.compartments = compartmentData
