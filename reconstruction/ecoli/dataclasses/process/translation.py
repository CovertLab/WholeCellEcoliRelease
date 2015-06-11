"""
SimulationData for translation process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/09/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
import numpy as np

class Translation(object):
	""" Translation """

	def __init__(self, raw_data, sim_data):
		self._buildMonomerData(raw_data, sim_data)
		self._buildTranslation(raw_data, sim_data)

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

			for aa in sim_data.amino_acid_1_to_3_ordered.viewkeys():
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
		fastRate = (np.log(2) / (2*units.min)).asUnit(1 / units.s)
		slowRate = (np.log(2) / (10*60*units.min)).asUnit(1 / units.s)

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

		degRate = np.zeros(len(raw_data.proteins))
		for i,m in enumerate(raw_data.proteins):
			if m['id'] not in ribosomalProteins:
				degRate[i] = NruleDegRate[m['seq'][0]].asNumber()
			else:
				degRate[i] = slowRate.asNumber()

		monomerData = np.zeros(
			size,
			dtype = [
				('id', 'a50'),
				('rnaId', 'a50'),
				('degRate', 'f8'),
				('length', 'i8'),
				('aaCounts', '{}i8'.format(nAAs)),
				('mw', 'f8'),
				('sequence', 'a{}'.format(maxSequenceLength)),
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
			'degRate'	:	1 / units.s,
			'length'	:	units.aa,
			'aaCounts'	:	units.aa,
			'mw'		:	units.g / units.mol,
			'sequence'  :   None
			}

		self.monomerData = UnitStructArray(monomerData, field_units)

	def _buildTranslation(self, raw_data, sim_data):
		from wholecell.utils.polymerize import PAD_VALUE

		sequences = self.monomerData["sequence"] # TODO: consider removing sequences

		maxLen = np.int64(
			self.monomerData["length"].asNumber().max()
			+ sim_data.constants.ribosomeElongationRate.asNumber(units.aa / units.s)
			)

		self.translationSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.translationSequences.fill(PAD_VALUE)

		aaIDs_singleLetter = sim_data.amino_acid_1_to_3_ordered.keys()

		aaMapping = {aa:i for i, aa in enumerate(aaIDs_singleLetter)}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.translationSequences[i, j] = aaMapping[letter]

		aaIDs = sim_data.amino_acid_1_to_3_ordered.values()

		self.translationMonomerWeights = (
			(
				sim_data.getter.getMass(aaIDs)
				- sim_data.getter.getMass(["WATER[c]"])
				)
			/ raw_data.constants['nAvogadro']
			).asNumber(units.fg)

		self.translationEndWeight = (sim_data.getter.getMass(["WATER[c]"]) / raw_data.constants['nAvogadro']).asNumber(units.fg)
