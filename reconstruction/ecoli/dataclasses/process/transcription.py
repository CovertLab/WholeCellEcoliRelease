"""
SimulationData for transcription process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import division

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
import numpy as np

class Transcription(object):
	""" Transcription """

	def __init__(self, raw_data, sim_data):
		self._buildRnaData(raw_data, sim_data)

	def _buildRnaData(self, raw_data, sim_data):
		assert all([len(rna['location']) == 1 for rna in raw_data.rnas])
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location'][0]) for rna in raw_data.rnas if len(rna['location']) == 1]
		rnaDegRates = np.log(2) / np.array([rna['halfLife'] for rna in raw_data.rnas]) # TODO: units
		rnaLens = np.array([len(rna['seq']) for rna in raw_data.rnas])
		ntCounts = np.array([
			(rna['seq'].count('A'), rna['seq'].count('C'),
				rna['seq'].count('G'), rna['seq'].count('U'))
			for rna in raw_data.rnas
			])
		expression = np.array([rna['expression'] for rna in raw_data.rnas])
		synthProb = expression * (
			np.log(2) / raw_data.parameters['cellCycleLen'].asNumber(units.s)
			+ rnaDegRates
			)
		
		synthProb /= synthProb.sum()

		mws = np.array([rna['mw'] for rna in raw_data.rnas]).sum(axis = 1)

		geneIds = np.array([rna['geneId'] for rna in raw_data.rnas])

		size = len(rnaIds)

		is23S = np.zeros(size, dtype = np.bool)
		is16S = np.zeros(size, dtype = np.bool)
		is5S = np.zeros(size, dtype = np.bool)

		for rnaIndex, rna in enumerate(raw_data.rnas):
			if rna["type"] == "rRNA" and rna["id"].startswith("RRL"):
				is23S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRS"):
				is16S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRF"):
				is5S[rnaIndex] = True

		sequences = [rna['seq'] for rna in raw_data.rnas]

		maxSequenceLength = max(len(sequence) for sequence in sequences)

		# TODO: Add units
		rnaData = np.zeros(
			size,
			dtype = [
				('id', 'a50'),
				# TODO: add expression to this table
				('synthProb', 'f8'),
				('degRate', 'f8'),
				('length', 'i8'),
				('countsACGU', '4i8'),
				('mw', 'f8'),
				('isMRna', 'bool'),
				('isMiscRna', 'bool'),
				('isRRna', 'bool'),
				('isTRna', 'bool'),
				('isRRna23S', 'bool'),
				('isRRna16S', 'bool'),
				('isRRna5S', 'bool'),
				('sequence', 'a{}'.format(maxSequenceLength)),
				('geneId', 'a50')
				]
			)

		rnaData['id'] = rnaIds
		rnaData['synthProb'] = synthProb
		rnaData['degRate'] = rnaDegRates
		rnaData['length'] = rnaLens
		rnaData['countsACGU'] = ntCounts
		rnaData['mw'] = mws
		rnaData['isMRna'] = [rna["type"] == "mRNA" for rna in raw_data.rnas]
		rnaData['isMiscRna'] = [rna["type"] == "miscRNA" for rna in raw_data.rnas]
		rnaData['isRRna'] = [rna["type"] == "rRNA" for rna in raw_data.rnas]
		rnaData['isTRna'] = [rna["type"] == "tRNA" for rna in raw_data.rnas]
		rnaData['isRRna23S'] = is23S
		rnaData['isRRna16S'] = is16S
		rnaData['isRRna5S'] = is5S
		rnaData['sequence'] = sequences
		rnaData['geneId'] = geneIds

		field_units = {
			'id'		:	None,
			'synthProb' :	None,
			'degRate'	:	1 / units.s,
			'length'	:	units.nt,
			'countsACGU':	units.nt,
			'mw'		:	units.g / units.mol,
			'isMRna'	:	None,
			'isMiscRna'	:	None,
			'isRRna'	:	None,
			'isTRna'	:	None,
			'isRRna23S'	:	None,
			'isRRna16S'	:	None,
			'isRRna5S'	:	None,
			'sequence'  :   None,
			'geneId'	:	None,
			}


		self.rnaData = UnitStructArray(rnaData, field_units)
		#self.getTrnaAbundanceData = getTrnaAbundanceAtGrowthRate