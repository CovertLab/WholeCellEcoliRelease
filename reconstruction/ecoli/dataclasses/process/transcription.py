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

#RNA_SEQ_ANALYSIS = "seal_rpkm"
RNA_SEQ_ANALYSIS = "rsem_tpm"

class Transcription(object):
	""" Transcription """

	def __init__(self, raw_data, sim_data):
		self._buildRnaData(raw_data, sim_data)
		self._buildTranscription(raw_data, sim_data)

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

		# Load expression from RNA-seq data
		expression = []
		for rna in raw_data.rnas:
			arb_exp = [x[sim_data.media_conditions] for x in eval("raw_data.rnaseq_{}_mean".format(RNA_SEQ_ANALYSIS)) if x['Gene'] == rna['geneId']]
			if len(arb_exp):
				expression.append(arb_exp[0])
			elif rna['type'] == 'mRNA' or rna['type'] == 'miscRNA':
				import ipdb; ipdb.set_trace()
				raise Exception('No RNA-seq data found for {}'.format(rna['id']))
			elif rna['type'] == 'rRNA' or rna['type'] == 'tRNA':
				expression.append(0.)
			else:
				raise Exception('Unknonw RNA {}'.format(rna['id']))

		expression = np.array(expression)
		synthProb = expression * (
			np.log(2) / sim_data.doubling_time.asNumber(units.s)
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
				('expression', 'float64'),
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
		rnaData["synthProb"] = synthProb
		rnaData["expression"] = expression
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
			'expression':	None,
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

	def _buildTranscription(self, raw_data, sim_data):
		from wholecell.utils.polymerize import PAD_VALUE

		sequences = self.rnaData["sequence"] # TODO: consider removing sequences

		maxLen = np.int64(
			self.rnaData["length"].asNumber().max()
			+ sim_data.constants.rnaPolymeraseElongationRate.asNumber(units.nt / units.s)
			)

		self.transcriptionSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.transcriptionSequences.fill(PAD_VALUE)

		ntMapping = {ntpId:i for i, ntpId in enumerate(["A", "C", "G", "U"])}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.transcriptionSequences[i, j] = ntMapping[letter]

		self.transcriptionMonomerWeights = (
			(
				sim_data.getter.getMass(sim_data.moleculeGroups.ntpIds)
				- sim_data.getter.getMass(["PPI[c]"])
				)
			/ raw_data.constants['nAvogadro']
			).asNumber(units.fg)

		self.transcriptionEndWeight = (sim_data.getter.getMass(["PPI[c]"]) / raw_data.constants['nAvogadro']).asNumber(units.fg)