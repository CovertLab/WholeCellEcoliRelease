"""
SimulationData for transcription process

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 03/06/2015
"""

from __future__ import division

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray
from wholecell.utils.polymerize import polymerize

RNA_SEQ_ANALYSIS = "rsem_tpm"
KCAT_ENDO_RNASE = 0.001
ESTIMATE_ENDO_RNASES = 5000
MAX_TIMESTEP_LEN = 2  # Determines length of padding values to add to transcript sequence matrix

class Transcription(object):
	"""
	SimulationData for the transcription process
	"""

	def __init__(self, raw_data, sim_data):
		self._buildRnaData(raw_data, sim_data)
		self._buildTranscription(raw_data, sim_data)

	def _buildRnaData(self, raw_data, sim_data):
		"""
		Build RNA-associated simulation data from raw data.
		"""

		assert all([len(rna['location']) == 1 for rna in raw_data.rnas])

		# Loads RNA IDs, degradation rates, lengths, and nucleotide compositions
		rnaIds = ['{}[{}]'.format(rna['id'], rna['location'][0])
            for rna in raw_data.rnas]
		rnaDegRates = np.log(2) / np.array([rna['halfLife'] for rna in raw_data.rnas]) # TODO: units
		rnaLens = np.array([len(rna['seq']) for rna in raw_data.rnas])

		ntCounts = np.array([
			(rna['seq'].count('A'), rna['seq'].count('C'),
			rna['seq'].count('G'), rna['seq'].count('U'))
			for rna in raw_data.rnas
			])

		# Load RNA expression from RNA-seq data
		expression = []

		for rna in raw_data.rnas:
			arb_exp = [x[sim_data.basal_expression_condition]
                for x in eval("raw_data.rna_seq_data.rnaseq_{}_mean".format(RNA_SEQ_ANALYSIS))
                if x['Gene'] == rna['geneId']]

			# If sequencing data is not found for rRNA or tRNA, initialize
            # expression to zero. For other RNA types, raise exception.
			if len(arb_exp) > 0:
				expression.append(arb_exp[0])
			elif rna['type'] == 'mRNA' or rna['type'] == 'miscRNA':
				raise Exception('No RNA-seq data found for {}'.format(rna['id']))
			elif rna['type'] == 'rRNA' or rna['type'] == 'tRNA':
				expression.append(0.)
			else:
				raise Exception('Unknown RNA {}'.format(rna['id']))

		expression = np.array(expression)

		# Calculate synthesis probabilities from expression and normalize
		synthProb = expression*(
			np.log(2) / sim_data.doubling_time.asNumber(units.s)
			+ rnaDegRates
			)
		synthProb /= synthProb.sum()

		# Calculate EndoRNase Km values
		Km = (KCAT_ENDO_RNASE*ESTIMATE_ENDO_RNASES/rnaDegRates) - expression

		# Load molecular weights and gene IDs
		mws = np.array([rna['mw'] for rna in raw_data.rnas]).sum(axis = 1)
		geneIds = np.array([rna['geneId'] for rna in raw_data.rnas])

		# Construct boolean arrays for rRNA types
		n_rnas = len(rnaIds)
		is23S = np.zeros(n_rnas, dtype = np.bool)
		is16S = np.zeros(n_rnas, dtype = np.bool)
		is5S = np.zeros(n_rnas, dtype = np.bool)

		for rnaIndex, rna in enumerate(raw_data.rnas):
			if rna["type"] == "rRNA" and rna["id"].startswith("RRL"):
				is23S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRS"):
				is16S[rnaIndex] = True

			if rna["type"] == "rRNA" and rna["id"].startswith("RRF"):
				is5S[rnaIndex] = True

		# Load sequence data
		sequences = [rna['seq'] for rna in raw_data.rnas]
		maxSequenceLength = max(len(sequence) for sequence in sequences)
		
		# Load IDs of protein monomers
		monomerIds = [rna['monomerId'] for rna in raw_data.rnas]

		oric_coordinate = raw_data.parameters['oriCCenter'].asNumber()
		terc_coordinate = raw_data.parameters['terCCenter'].asNumber()
		genome_length = len(raw_data.genome_sequence)

		def get_replication_coordinate(coordinate):
			replication_coordinate = ((coordinate - terc_coordinate)
				% genome_length + terc_coordinate - oric_coordinate
				)

			if replication_coordinate < 0:
				replication_coordinate += 1

			return replication_coordinate

		replicationCoordinate = [
			get_replication_coordinate(x["coordinate"])
			for x in raw_data.rnas]

		rnaData = np.zeros(
			n_rnas,
			dtype = [
				('id', 'a50'),
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
				('isRProtein', 'bool'),
				('isRnap',	'bool'),
				('sequence', 'a{}'.format(maxSequenceLength)),
				('geneId', 'a50'),
				('KmEndoRNase', 'f8'),
				('replicationCoordinate', 'int64'),
				]
			)

		rnaData['id'] = rnaIds
		rnaData['degRate'] = rnaDegRates
		rnaData['length'] = rnaLens
		rnaData['countsACGU'] = ntCounts
		rnaData['mw'] = mws
		rnaData['isMRna'] = [rna["type"] == "mRNA" for rna in raw_data.rnas]
		rnaData['isMiscRna'] = [rna["type"] == "miscRNA" for rna in raw_data.rnas]
		rnaData['isRRna'] = [rna["type"] == "rRNA" for rna in raw_data.rnas]
		rnaData['isTRna'] = [rna["type"] == "tRNA" for rna in raw_data.rnas]
		rnaData['isRProtein'] = [
            "{}[c]".format(x) in sim_data.moleculeGroups.rProteins
            for x in monomerIds]
		rnaData['isRnap'] = [
            "{}[c]".format(x) in sim_data.moleculeGroups.rnapIds
            for x in monomerIds]
		rnaData['isRRna23S'] = is23S
		rnaData['isRRna16S'] = is16S
		rnaData['isRRna5S'] = is5S
		rnaData['sequence'] = sequences
		rnaData['geneId'] = geneIds
		rnaData['KmEndoRNase'] = Km
		rnaData['replicationCoordinate'] = replicationCoordinate

		field_units = {
			'id': None,
			'degRate': 1 / units.s,
			'length': units.nt,
			'countsACGU': units.nt,
			'mw': units.g / units.mol,
			'isMRna': None,
			'isMiscRna': None,
			'isRRna': None,
			'isTRna': None,
			'isRRna23S': None,
			'isRRna16S': None,
			'isRRna5S':	None,
			'isRProtein': None,
			'isRnap': None,
			'sequence': None,
			'geneId': None,
			'KmEndoRNase': units.mol / units.L,
			'replicationCoordinate': None,
			}

		self.rnaExpression = {}
		self.rnaSynthProb = {}

		# Set basal expression and synthesis probabilities - conditional values
        # are set in the fitter.
		self.rnaExpression["basal"] = expression / expression.sum()
		self.rnaSynthProb["basal"] = synthProb / synthProb.sum()

		self.rnaData = UnitStructArray(rnaData, field_units)


	def _buildTranscription(self, raw_data, sim_data):
		"""
		Build transcription-associated simulation data from raw data.
		"""
		sequences = self.rnaData["sequence"] # TODO: consider removing sequences

		# Construct transcription sequence matrix
		maxLen = np.int64(
			self.rnaData["length"].asNumber().max()
			+ MAX_TIMESTEP_LEN*sim_data.growthRateParameters.rnaPolymeraseElongationRate.asNumber(units.nt/units.s)
			)

		self.transcriptionSequences = np.empty((sequences.shape[0], maxLen), np.int8)
		self.transcriptionSequences.fill(polymerize.PAD_VALUE)

		ntMapping = {ntpId: i for i, ntpId in enumerate(["A", "C", "G", "U"])}

		for i, sequence in enumerate(sequences):
			for j, letter in enumerate(sequence):
				self.transcriptionSequences[i, j] = ntMapping[letter]

		# Calculate weights of transcript nucleotide monomers
		self.transcriptionMonomerWeights = (
			(
				sim_data.getter.getMass(sim_data.moleculeGroups.ntpIds)
				- sim_data.getter.getMass(["PPI[c]"])
				)
			/ raw_data.constants['nAvogadro']
			).asNumber(units.fg)

		self.transcriptionEndWeight = ((sim_data.getter.getMass(["PPI[c]"])
            / raw_data.constants['nAvogadro']).asNumber(units.fg))
