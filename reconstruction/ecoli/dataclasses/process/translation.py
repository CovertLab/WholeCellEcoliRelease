"""
SimulationData for translation process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/09/2015
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import six

from wholecell.sim.simulation import MAX_TIME_STEP
from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.polymerize import polymerize
from wholecell.utils.random import make_elongation_rates


PROCESS_MAX_TIME_STEP = 2.


class Translation(object):
	""" Translation """

	def __init__(self, raw_data, sim_data):
		self.max_time_step = min(MAX_TIME_STEP, PROCESS_MAX_TIME_STEP)

		self._buildMonomerData(raw_data, sim_data)
		self._buildTranslation(raw_data, sim_data)
		self._buildTranslationEfficiency(raw_data, sim_data)
		self._build_elongation_rates(raw_data, sim_data)

	def _buildMonomerData(self, raw_data, sim_data):
		assert all([len(protein['location']) == 1 for protein in raw_data.proteins])
		ids = ['{}[{}]'.format(protein['id'], protein['location'][0]) for protein in raw_data.proteins]

		rnaIds = []

		for protein in raw_data.proteins:
			rnaId = protein['rnaId']

			rnaLocation = None
			for rna in raw_data.rnas:
				if rna['id'] == rnaId:
					assert len(rna['location']) == 1
					rnaLocation = rna['location'][0]
					break

			rnaIds.append('{}[{}]'.format(
				rnaId,
				rnaLocation
				))

		lengths = []
		aaCounts = []
		sequences = []

		for protein in raw_data.proteins:
			sequence = protein['seq']

			counts = []

			for aa in sim_data.amino_acid_1_to_3_ordered:
				counts.append(
					sequence.count(aa)
					)

			lengths.append(len(sequence))
			aaCounts.append(counts)
			sequences.append(sequence)

		maxSequenceLength = max(len(seq) for seq in sequences)

		mws = np.array([protein['mw'] for protein in raw_data.proteins]).sum(axis = 1)

		size = len(rnaIds)

		nAAs = len(aaCounts[0])

		# Calculate degradation rates based on N-rule
		# TODO: citation
		deg_rate_units = 1 / units.s
		fastRate = (np.log(2) / (2*units.min)).asNumber(deg_rate_units)
		slowRate = (np.log(2) / (10*60*units.min)).asNumber(deg_rate_units)

		fastAAs = ["R", "K", "F", "L", "W", "Y"]
		slowAAs = ["H", "I", "D", "E", "N", "Q", "C", "A", "S", "T", "G", "V", "M"]
		noDataAAs = ["P", "U"]

		NruleDegRate = {}
		NruleDegRate.update(
			(fastAA, fastRate) for fastAA in fastAAs
			)
		NruleDegRate.update(
			(slowAA, slowRate) for slowAA in slowAAs
			)
		NruleDegRate.update(
			(noDataAA, slowRate) for noDataAA in noDataAAs
			) # Assumed slow rate because of no data

		# Build list of ribosomal proteins
		# Give all ribosomal proteins the slowAA rule
		ribosomalProteins = []
		ribosomalProteins.extend([x[:-3] for x in sim_data.moleculeGroups.s30_proteins])
		ribosomalProteins.extend([x[:-3] for x in sim_data.moleculeGroups.s50_proteins])

		# Get degradation rates from measured protein half lives
		measured_deg_rates = {
			p['id']: (np.log(2) / p['half life']).asNumber(deg_rate_units)
			for p in raw_data.protein_half_lives
			}

		degRate = np.zeros(len(raw_data.proteins))
		for i,m in enumerate(raw_data.proteins):
			if m['id'] in measured_deg_rates:
				degRate[i] = measured_deg_rates[m['id']]
			elif m['id'] not in ribosomalProteins:
				degRate[i] = NruleDegRate[m['seq'][0]]
			else:
				degRate[i] = slowRate

		monomerData = np.zeros(
			size,
			dtype = [
				('id', 'U50'),
				('rnaId', 'U50'),
				('degRate', 'f8'),
				('length', 'i8'),
				('aaCounts', '{}i8'.format(nAAs)),
				('mw', 'f8'),
				('sequence', 'U{}'.format(maxSequenceLength)),
				]
			)

		monomerData['id'] = ids
		monomerData['rnaId'] = rnaIds
		monomerData['degRate'] = degRate
		monomerData['length'] = lengths
		monomerData['aaCounts'] = aaCounts
		monomerData['mw'] = mws
		monomerData['sequence'] = sequences

		field_units = {
			'id'		:	None,
			'rnaId'		:	None,
			'degRate'	:	deg_rate_units,
			'length'	:	units.aa,
			'aaCounts'	:	units.aa,
			'mw'		:	units.g / units.mol,
			'sequence'  :   None
			}

		self.monomerData = UnitStructArray(monomerData, field_units)
		self.n_monomers = len(self.monomerData)

	def _buildTranslation(self, raw_data, sim_data):
		sequences = self.monomerData["sequence"] # TODO: consider removing sequences

		maxLen = np.int64(
			self.monomerData["length"].asNumber().max()
			+ self.max_time_step * sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s)
			)

		self.translationSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.translationSequences.fill(polymerize.PAD_VALUE)

		aaIDs_singleLetter = six.viewkeys(sim_data.amino_acid_1_to_3_ordered)

		aaMapping = {aa:i for i, aa in enumerate(aaIDs_singleLetter)}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.translationSequences[i, j] = aaMapping[letter]

		aaIDs = list(sim_data.amino_acid_1_to_3_ordered.values())

		self.translationMonomerWeights = (
			(
				sim_data.getter.getMass(aaIDs)
				- sim_data.getter.getMass([sim_data.moleculeIds.water])
				)
			/ sim_data.constants.nAvogadro
			).asNumber(units.fg)

		self.translationEndWeight = (sim_data.getter.getMass([sim_data.moleculeIds.water]) / sim_data.constants.nAvogadro).asNumber(units.fg)

	def _buildTranslationEfficiency(self, raw_data, sim_data):
		monomerIds = [
			x["id"] + sim_data.getter.get_location_tag(x["id"])
			for x in raw_data.proteins]
		monomerIdToGeneId = {
			x["id"] + sim_data.getter.get_location_tag(x["id"]): x["geneId"]
			for x in raw_data.proteins}
		geneIdToTrEff = {
			x["geneId"]: x["translationEfficiency"]
			for x in raw_data.translationEfficiency
			if type(x["translationEfficiency"]) == float}
		trEffs = []
		for monomerId in monomerIds:
			geneId = monomerIdToGeneId[monomerId]
			if geneId in geneIdToTrEff:
				trEffs.append(geneIdToTrEff[geneId])
			else:
				trEffs.append(np.nan)

		self.translationEfficienciesByMonomer = np.array(trEffs)
		self.translationEfficienciesByMonomer[np.isnan(self.translationEfficienciesByMonomer)] = np.nanmean(self.translationEfficienciesByMonomer)


	def _build_elongation_rates(self, raw_data, sim_data):
		protein_ids = self.monomerData['id']
		ribosomal_protein_ids = sim_data.moleculeGroups.rProteins

		protein_indexes = {
			protein: index
			for index, protein in enumerate(protein_ids)}

		ribosomal_proteins = {
			rprotein: protein_indexes.get(rprotein, -1)
			for rprotein in ribosomal_protein_ids}

		self.rprotein_indexes = np.array([
			index
			for index in ribosomal_proteins.values()
			if index >= 0], dtype=np.int64)

		self.basal_elongation_rate = sim_data.constants.ribosomeElongationRateBasal.asNumber(units.aa / units.s)
		self.max_elongation_rate = sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s)
		self.elongation_rates = np.full(
			self.n_monomers,
			self.basal_elongation_rate,
			dtype=np.int64)

		self.elongation_rates[self.rprotein_indexes] = self.max_elongation_rate

	def make_elongation_rates(
			self,
			random,
			base,
			time_step,
			variable_elongation=False):

		return make_elongation_rates(
			random,
			self.n_monomers,
			base,
			self.rprotein_indexes,
			self.max_elongation_rate,
			time_step,
			variable_elongation)
