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
		# The attributes of active RNA polymerases are given as:
		# - domain_index (32-bit int): Domain index of the chromosome domain
		# that the RNAP is bound to. This value is used to split the RNAPs at
		# cell division.
		# - coordinates (64-bit int): Location of the RNAP on the chromosome,
		# in base pairs from origin.
		# - direction (bool): True if RNAP is moving in the positive direction
		# of the coordinate, False if RNAP is moving in the negative direction.
		# This is determined by the orientation of the gene that the RNAP is
		# transcribing.
		RNAP_mass = self.bulkMolecules.bulkData['mass'][
			self.bulkMolecules.bulkData['id'] == sim_data.moleculeIds.rnapFull]
		RNAP_attributes = {
			'domain_index': 'i4',
			'coordinates': 'i8',
			'direction': '?',
			}

		self.uniqueMolecules.addToUniqueState('active_RNAP', RNAP_attributes, RNAP_mass)

		# RNAPs are divided based on the index of the chromosome domain they
		# are bound to
		sim_data.moleculeGroups.unique_molecules_domain_index_division.append(
			'active_RNAP')

		# Add RNAs
		# This molecule represents all RNAs that should be represented as
		# unique molecules. These RNAs include all partially transcribed RNAs
		# and fully transcribed mRNAs. The attributes of RNAs are given as:
		# - TU_index (64-bit int): Index of the transcription unit that the
		# RNA molecule is representing. Determines the sequence and the length
		# of the fully elongated RNA.
		# - transcript_length (64-bit int): Current length of the RNA.
		# - is_mRNA (bool): True if RNA represents an mRNA molecule.
		# - is_full_transcript (bool): True if RNA is fully transcribed and
		# released from the RNA polymerase, False if RNA is being actively
		# transcribed. This cannot be True if 'is_mRNA' is False, since all
		# non-mRNAs are represented as bulk molecules when they are fully
		# transcribed.
		# - can_translate (bool): True if the 5' end of the mRNA molecule is
		# available for translation. This flag is used only for mRNA molecules.
		# - RNAP_index (64-bit int): Unique index of the RNA polymerase that is
		# synthesizing the RNA or have synthesized the mRNA. For fully
		# transcribed mRNAs that are added at initialization, this attribute
		# is set to -1.
		RNA_mass = (units.g/units.mol) * np.zeros_like(RNAP_mass)
		RNA_attributes = {
			'TU_index': 'i8',
			'transcript_length': 'i8',
			'is_mRNA': '?',
			'is_full_transcript': '?',
			'can_translate': '?',
			'RNAP_index': 'i8',
			}

		self.uniqueMolecules.addToUniqueState('RNA', RNA_attributes, RNA_mass)

		# Fully transcribed mRNAs are divided binomially, partial transcripts
		# are divided based on which chromosome domains their associated RNAPs
		# are bound to
		sim_data.moleculeGroups.unique_molecules_RNA_division.append('RNA')

		# Add active ribosomes
		# The attributes of active ribosomes are given as:
		# - protein_index (64-bit int): Index of the protein monomer that the
		# ribosome is translating
		# - peptide_length (64-bit int): Current length of the polypeptide that
		# the ribosome is translating, in number of amino acids
		# - mRNA_index (64-bit int): Unique index of the mRNA that the ribosome
		# is bound to
		# - pos_on_mRNA (64-bit int): Location of the ribosome on the bound
		# mRNA, in number of bases from the transcription start site
		# TODO: This is a bad hack that works because in the parca
		# I have forced expression to be these subunits only
		ribosome_30S_mass = self.bulkMolecules.bulkData['mass'][
			self.bulkMolecules.bulkData['id'] == sim_data.moleculeIds.s30_fullComplex]
		ribosome_50S_mass = self.bulkMolecules.bulkData['mass'][
			self.bulkMolecules.bulkData['id'] == sim_data.moleculeIds.s50_fullComplex]
		ribosome_mass = ribosome_30S_mass + ribosome_50S_mass
		ribosome_attributes = {
			'protein_index': 'i8',
			'peptide_length': 'i8',
			'mRNA_index': 'i8',
			'pos_on_mRNA': 'i8',
			}
		self.uniqueMolecules.addToUniqueState('active_ribosome', ribosome_attributes, ribosome_mass)

		# Active ribosomes are divided such that they always follow the mRNA
		# molecule they are bound to
		sim_data.moleculeGroups.unique_molecules_active_ribosome_division.append(
			'active_ribosome')

		# Add full chromosomes
		# One full chromosome molecule is added when chromosome replication is
		# complete, and sets cell division to happen after a length of time
		# specified by the D period (if D_PERIOD_DIVISION is set to True).
		# The 'has_triggered_division' attribute is initially set to False, and
		# is reset to True when division_time was reached and the cell has
		# divided. The 'domain_index' keeps track of the index of the oldest
		# chromosome domain that is part of the full chromosome.
		fullChromosomeMass = (units.g/units.mol) * (
			stateFunctions.createMassesByCompartments(raw_data.full_chromosome))
		fullChromosomeAttributes = {
			'division_time': 'f8',
			'has_triggered_division': '?',
			'domain_index': 'i4',
			}

		self.uniqueMolecules.addToUniqueState('full_chromosome', fullChromosomeAttributes, fullChromosomeMass)

		# Full chromosomes are divided based on their domain index
		sim_data.moleculeGroups.unique_molecules_domain_index_division.append(
			'full_chromosome')


		# Add chromosome domains
		# Chromosome domains are zero-mass molecules that accounts for the
		# structures of replicating chromosomes. Each replication initiation
		# event creates two new chromosome domains that are given a unique
		# integer 'domain_index'. These two new domains are child domains of
		# the original domain that the origin belonged to.
		chromosome_domain_mass = (units.g/units.mol) * np.zeros_like(RNAP_mass)
		chromosome_domain_attributes = {
			'domain_index': 'i4',
			'child_domains': ('i4', 2)
			}

		self.uniqueMolecules.addToUniqueState('chromosome_domain', chromosome_domain_attributes, chromosome_domain_mass)

		# Chromosome domains are divided based on their domain index
		sim_data.moleculeGroups.unique_molecules_domain_index_division.append(
			'chromosome_domain')


		# Add active replisomes
		# Note that the replisome does not functionally replicate the
		# chromosome, but instead keeps track of the mass associated with
		# essential subunits of the replisome complex. The list of essential
		# subunits and their stoichiometry were taken from Reyes-Lamothe et
		# al., 2010.
		trimer_ids = sim_data.moleculeGroups.replisome_trimer_subunits
		monomer_ids = sim_data.moleculeGroups.replisome_monomer_subunits

		trimer_mass = [self.bulkMolecules.bulkData['mass'][
			self.bulkMolecules.bulkData['id'] == id].asNumber(units.g/units.mol)
			for id in trimer_ids]
		monomer_mass = [self.bulkMolecules.bulkData['mass'][
			self.bulkMolecules.bulkData['id'] == id].asNumber(units.g/units.mol)
			for id in monomer_ids]

		replisomeMass = (units.g/units.mol) * (
				3*np.sum(trimer_mass, axis=0) + np.sum(monomer_mass, axis=0))

		replisomeAttributes = {
			'domain_index': 'i4',
			'right_replichore': '?',
			'coordinates': 'i8',
			}

		self.uniqueMolecules.addToUniqueState('active_replisome', replisomeAttributes, replisomeMass)

		# Active replisomes are divided based on their domain index
		sim_data.moleculeGroups.unique_molecules_domain_index_division.append(
			'active_replisome')


		# Add origins of replication
		# Note that origins are conceptual molecules and have zero mass. The
		# chromosomeIndexes of oriC's determine the chromosomeIndexes of the
		# new partial chromosomes and replisomes initiated on the same oriC.
		originMass = (units.g/units.mol) * np.zeros_like(RNAP_mass)
		originAttributes = {
			'domain_index': 'i4',
			}

		self.uniqueMolecules.addToUniqueState('oriC', originAttributes, originMass)

		# oriC's are divided based on their domain index
		sim_data.moleculeGroups.unique_molecules_domain_index_division.append(
			'oriC')

		# Add promoters
		# Promoters are sequences on the DNA where RNA polymerases bind to and
		# initiate transcription. They can also be bound to transcription
		# factors(TFs), if the transcription unit associated with the promoter
		# is regulated by TFs. The promoter itself has zero mass but can hold
		# the mass of the transcription factor that it is bound to. Its
		# attributes are given as:
		# - TU_index (64-bit int): Index of the transcription unit that
		# the promoter is associated with. This determines which TFs the
		# promoter can bind to, and how the transcription probability is
		# calculated.
		# - coordinates (64-bit int): Location of the promoter on the
		# chromosome, in base pairs from origin. This value does not change
		# after the molecule is initialized.
		# - domain_index (32-bit int): Domain index of the chromosome domain
		# that the promoter belongs to. This value is used to split the
		# promoters at cell division.
		# - bound_TF (boolean array of length n_tf): A boolean array that
		# shows which TFs the promoter is bound to. Note that one promoter can
		# bind to multiple TFs.
		n_tf = len(sim_data.process.transcription_regulation.tf_ids)

		promoter_mass = (units.g/units.mol) * np.zeros_like(RNAP_mass)
		promoter_attributes = {
			'TU_index': 'i8',
			'coordinates': 'i8',
			'domain_index': 'i4',
			'bound_TF': ('?', n_tf),
			}

		self.uniqueMolecules.addToUniqueState('promoter', promoter_attributes, promoter_mass)

		# Promoters are divided based on their domain index
		sim_data.moleculeGroups.unique_molecules_domain_index_division.append(
			'promoter'
			)

		# Add DnaA boxes
		# DnaA boxes are 9-base sequence motifs on the DNA that bind to the
		# protein DnaA. Except for DnaA boxes close to the origin, these boxes
		# serve no functional role in replication initiation, but can
		# effectively titrate away free DnaA molecules and control its
		# concentration. The molecule itself has zero mass but it can hold the
		# mass of the DnaA protein that it is bound to. Its attributes are
		# given as:
		# - coordinates (64-bit int): Location of the middle base (5th base) of
		# the DnaA box, in base pairs from origin. This value does not change
		# after the molecule is initialized.
		# - domain_index (32-bit int): Domain index of the chromosome domain
		# that the DnaA box belongs to. This value is used to allocate DnaA
		# boxes to the two daughter cells at cell division.
		# - DnaA_bound (boolean): True if bound to a DnaA protein, False if not
		DnaA_box_mass = (units.g/units.mol) * np.zeros_like(RNAP_mass)
		DnaA_box_attributes = {
			'coordinates': 'i8',
			'domain_index': 'i4',
			'DnaA_bound': '?',
			}

		self.uniqueMolecules.addToUniqueState('DnaA_box',
			DnaA_box_attributes, DnaA_box_mass)

		# DnaA boxes are divided based on their domain index
		sim_data.moleculeGroups.unique_molecules_domain_index_division.append(
			'DnaA_box')

	def _buildCompartments(self, raw_data, sim_data):
		compartmentData = np.empty(len(raw_data.compartments),
			dtype = [('id','a20'),('compartmentAbbreviation', 'a1')])

		compartmentData['id'] = [x['id'] for x in raw_data.compartments]
		compartmentData['compartmentAbbreviation'] = [x['abbrev'] for x in raw_data.compartments]
		self.compartments = compartmentData
