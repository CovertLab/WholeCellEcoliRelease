"""
SimulationData moleculeIds
"""

from __future__ import absolute_import, division, print_function

class MoleculeIds(object):
	"""
	Helper class to extract molecule IDs of "special" molecules. All values
	returned are strings.
	"""

	def __init__(self, raw_data, sim_data):
		self._buildMoleculeIds()

	def _buildMoleculeIds(self):
		molecule_ids = {
			'full_RNAP': 'APORNAP-CPLX[c]',
			's30_full_complex': 'CPLX0-3953[c]',
			's50_full_complex': 'CPLX0-3962[c]',
			'DnaA': 'PD03831[c]',
			'DnaA_ATP_complex': 'MONOMER0-160[c]',
			'LPS': 'CPD0-939[c]',
			'murein': 'CPD-12261[p]',
			'glycogen': 'glycogen-monomer[c]',
			'ppGpp': 'GUANOSINE-5DP-3DP[c]',
			'oriC_site': 'G0-10506',
			'terC_site': 'TERC',
			'RelA': 'RELA-MONOMER[c]',
			'SpoT': 'SPOT-MONOMER[c]',
			'water': 'WATER[c]',
			'proton': 'PROTON[c]',
			'ppi': 'PPI[c]',
			'full_chromosome': 'CHROM_FULL[c]'
		}

		self.__dict__.update(molecule_ids)
