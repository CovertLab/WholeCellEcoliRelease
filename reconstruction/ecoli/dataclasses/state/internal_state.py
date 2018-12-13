"""
SimulationData state associated data

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

from reconstruction.ecoli.dataclasses.state.bulkMolecules import BulkMolecules
from reconstruction.ecoli.dataclasses.state.uniqueMolecules import UniqueMolecules
from reconstruction.ecoli.dataclasses.state import stateFunctions

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
		"""
		Add data (IDs and mass) for all classes of bulk molecules.
		"""

		# Set metabolites
		metaboliteIds = stateFunctions.createIdsWithCompartments(raw_data.metabolites)
		metaboliteMasses = (units.g/units.mol) * (
			stateFunctions.createMetaboliteMassesByCompartments(raw_data.metabolites, 7, 11))

		self.bulkMolecules.addToBulkState(metaboliteIds, metaboliteMasses)

		# Set water
		waterIds = stateFunctions.createIdsWithCompartments(raw_data.water)
		waterMasses = (units.g/units.mol) * (
			stateFunctions.createMetaboliteMassesByCompartments(raw_data.water, 8, 11))

		self.bulkMolecules.addToBulkState(waterIds, waterMasses)

		# Set RNA
		rnaIds = stateFunctions.createIdsWithCompartments(raw_data.rnas)
		rnaMasses = (units.g/units.mol) * (
			stateFunctions.createMassesByCompartments(raw_data.rnas))

		self.bulkMolecules.addToBulkState(rnaIds, rnaMasses)

		# Set proteins
		proteinIds = stateFunctions.createIdsWithCompartments(raw_data.proteins)
		proteinMasses = (units.g/units.mol) * (
			stateFunctions.createMassesByCompartments(raw_data.proteins))

		self.bulkMolecules.addToBulkState(proteinIds, proteinMasses)

		# Set complexes
		complexIds = stateFunctions.createIdsWithCompartments(raw_data.proteinComplexes)
		complexMasses = (units.g/units.mol) * (
			stateFunctions.createMassesByCompartments(raw_data.proteinComplexes))

		self.bulkMolecules.addToBulkState(complexIds, complexMasses)

		# Set modified forms
		modifiedFormIds = stateFunctions.createIdsWithCompartments(raw_data.modifiedForms)
		modifiedFormMasses = (units.g/units.mol) * (
			stateFunctions.createModifiedFormMassesByCompartments(raw_data.modifiedForms))

		self.bulkMolecules.addToBulkState(modifiedFormIds, modifiedFormMasses)

		# Set fragments
		fragments = []
		
		for x in raw_data.polymerized:
			if x['is_ntp']:
				if not x['is_end']:
					temp = x
					temp['id'] = x['id'].replace('Polymerized','Fragment')
					fragments.append(temp)
					
		fragmentsIds = stateFunctions.createIdsWithCompartments(fragments)
		fragmentsMasses = (units.g/units.mol) * (
			stateFunctions.createMassesByCompartments(fragments))

		self.bulkMolecules.addToBulkState(fragmentsIds, fragmentsMasses)


	def _buildUniqueMolecules(self, raw_data, sim_data):
		"""
		Add data (name, mass, and attribute data structure) for all classes of
		unique molecules.
		"""

		# Add active RNA polymerase
		rnaPolyComplexMass = self.bulkMolecules.bulkData["mass"][
			self.bulkMolecules.bulkData["id"] == sim_data.moleculeIds.rnapFull]
		rnaPolyAttributes = {
			"rnaIndex": "i8",
			"transcriptLength": "i8"
			}

		self.uniqueMolecules.addToUniqueState('activeRnaPoly', rnaPolyAttributes, rnaPolyComplexMass)

		# Add active ribosome
		# TODO: This is a bad hack that works because in the fitter
		# I have forced expression to be these subunits only
		ribosome30SMass = self.bulkMolecules.bulkData["mass"][
			self.bulkMolecules.bulkData["id"] == sim_data.moleculeIds.s30_fullComplex]
		ribosome50SMass = self.bulkMolecules.bulkData["mass"][
			self.bulkMolecules.bulkData["id"] == sim_data.moleculeIds.s50_fullComplex]
		ribosomeMass = ribosome30SMass + ribosome50SMass
		ribosomeAttributes = {
			"proteinIndex": "i8",
			"peptideLength": "i8",
			}
		self.uniqueMolecules.addToUniqueState("activeRibosome", ribosomeAttributes, ribosomeMass)


		# Add full chromosomes
		# One full chromosome molecule is added when chromosome replication is
		# complete, and sets cell division to happen after a length of time
		# specified by the D period (if D_PERIOD_DIVISION is set to True).
		# The "has_induced_division" attribute is initially set to False, and
		# is reset to True when division_time was reached and the cell has
		# divided. The "mother_domain_index" keeps track of the index of the
		# oldest chromosome domain of the full chromosome.
		fullChromosomeMass = (units.g/units.mol) * (
			stateFunctions.createMassesByCompartments(raw_data.full_chromosome))
		fullChromosomeAttributes = {
			"division_time": "f8",
			"has_induced_division": "?",
			"mother_domain_index": "i8",
			}

		self.uniqueMolecules.addToUniqueState('fullChromosome', fullChromosomeAttributes, fullChromosomeMass)


		# Add chromosome domains
		# Chromosome domains are zero-mass molecules that accounts for the
		# structures of replicating chromosomes. Each replication initiation
		# event creates two new chromosome domains that are given a unique
		# integer "domain_index". These two new domains are child domains of
		# the original domain that the origin belonged to.
		chromosome_domain_mass = (units.g/units.mol) * np.zeros_like(rnaPolyComplexMass)
		chromosome_domain_attributes = {
			"domain_index": "i8",
			"child_domains": ("i8", 2)
			}

		# Placeholder value for domains without children domains
		sim_data.process.replication.no_child_place_holder = -1

		self.uniqueMolecules.addToUniqueState('chromosome_domain', chromosome_domain_attributes, chromosome_domain_mass)


		# Add active replisomes
		# Note that the replisome does not functionally replicate the
		# chromosome, but instead keeps track of the mass associated with
		# essential subunits of the replisome complex. The list of essential
		# subunits and their stoichiometry were taken from Reyes-Lamothe et
		# al., 2010.
		trimer_ids = sim_data.moleculeGroups.replisome_trimer_subunits
		monomer_ids = sim_data.moleculeGroups.replisome_monomer_subunits

		trimer_mass = [self.bulkMolecules.bulkData["mass"][
			self.bulkMolecules.bulkData["id"] == id].asNumber(units.g/units.mol)
			for id in trimer_ids]
		monomer_mass = [self.bulkMolecules.bulkData["mass"][
			self.bulkMolecules.bulkData["id"] == id].asNumber(units.g/units.mol)
			for id in monomer_ids]

		replisomeMass = (units.g/units.mol) * (
				3*np.sum(trimer_mass, axis=0) + np.sum(monomer_mass, axis=0))

		replisomeAttributes = {
			"domain_index": "i8",
			"right_replichore": "?",
			"coordinates": "i8",
			}

		self.uniqueMolecules.addToUniqueState('active_replisome', replisomeAttributes, replisomeMass)


		# Add origins of replication
		# Note that origins are conceptual molecules and have zero mass. The
		# chromosomeIndexes of oriC's determine the chromosomeIndexes of the
		# new partial chromosomes and replisomes initiated on the same oriC.
		originMass = (units.g/units.mol) * np.zeros_like(rnaPolyComplexMass)
		originAttributes = {
			"domain_index": "i8",
			}

		self.uniqueMolecules.addToUniqueState('originOfReplication', originAttributes, originMass)


	def _buildCompartments(self, raw_data, sim_data):
		compartmentData = np.empty(len(raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in raw_data.compartments]
		self.compartments = compartmentData
