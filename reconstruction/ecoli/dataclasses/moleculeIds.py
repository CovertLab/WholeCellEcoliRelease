"""
SimulationData moleculeIds

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 07/18/2018
"""

from __future__ import division

class MoleculeIds(object):
	"""
	Helper class to extract molecule IDs of "special" molecules. All values
	returned are strings.
	"""

	def __init__(self, raw_data, sim_data):
		self._buildMoleculeIds()

	def _buildMoleculeIds(self):
		moleculeIds = {
			'rnapFull':	'APORNAP-CPLX[c]',
			's30_fullComplex': 'CPLX0-3953[c]',
			's50_fullComplex': 'CPLX0-3962[c]',
			'fullChromosome': "CHROM_FULL[c]",
		}

		self.__dict__.update(moleculeIds)
