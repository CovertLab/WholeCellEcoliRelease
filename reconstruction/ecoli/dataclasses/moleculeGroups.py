"""
SimulationData molecule groups

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division


class moleculeGroups(object):
	""" moleculeGroups """

	def __init__(self, raw_data, sim_data):
		self._buildMoleculeGroups()

	def _buildMoleculeGroups(self):
		moleculeGroups = {
			'ntpIds'			:	["ATP[c]","CTP[c]","GTP[c]","UTP[c]"],
			'dNtpIds'			:	["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"],
			'dNmpIds'			:	["DAMP[c]", "DCMP[c]", "DGMP[c]", "DTMP[c]"],
			'dNmpNuclearIds'	:	["DAMP[n]", "DCMP[n]", "DGMP[n]", "DTMP[n]"],
			'rnapIds'			:	["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"],
			#'polymerizedAA_IDs'	:	self._polymerizedAA_IDs, # TODO: end weight # TODO: Add these groups!
			#'polymerizedNT_IDs'	:	self._polymerizedNT_IDs, # TODO: end weight # TODO: Add these groups!
			#'polymerizedDNT_IDs':	self._polymerizedDNT_IDs, # TODO: Add these groups!
			's30_proteins'			:	['EG10912-MONOMER[c]', 'EG10916-MONOMER[c]', 'EG10906-MONOMER[c]',
										'EG10914-MONOMER[c]', 'EG10909-MONOMER[c]', 'EG10903-MONOMER[c]',
										'EG10911-MONOMER[c]', 'EG10904-MONOMER[c]', 'EG10900-MONOMER[c]',
										'EG10901-MONOMER[c]', 'EG10905-MONOMER[c]',
										'EG10915-MONOMER[c]', 'EG10918-MONOMER[c]', 'EG10919-MONOMER[c]',
										'EG10907-MONOMER[c]', 'EG11508-MONOMER[c]', 'EG10908-MONOMER[c]',
										'EG10920-MONOMER[c]', 'EG10910-MONOMER[c]', 'EG10902-MONOMER[c]',
										'EG10917-MONOMER[c]', 'EG10913-MONOMER[c]'],
			's30_16sRRNA'			:	['RRSA-RRNA[c]','RRSB-RRNA[c]','RRSC-RRNA[c]','RRSD-RRNA[c]',
										'RRSE-RRNA[c]','RRSG-RRNA[c]','RRSH-RRNA[c]'],
			's30_fullComplex'		:	['CPLX-30SA[c]'],
			's50_proteins'			:	['EG10872-MONOMER[c]', 'EG10879-MONOMER[c]',
										'EG11232-MONOMER[c]', 'EG10877-MONOMER[c]',
										'EG10876-MONOMER[c]', 'EG10892-MONOMER[c]', 'EG10874-MONOMER[c]',
										'EG50001-MONOMER[c]', 'EG10875-MONOMER[c]', 'EG10884-MONOMER[c]',
										'EG11231-MONOMER[c]', 'EG10887-MONOMER[c]',
										'EG10878-MONOMER[c]', 'EG10886-MONOMER[c]', 'EG10870-MONOMER[c]',
										'EG10889-MONOMER[c]', 'EG10891-MONOMER[c]', 'EG10888-MONOMER[c]',
										'EG50002-MONOMER[c]', 'EG10869-MONOMER[c]', 'EG10882-MONOMER[c]',
										'EG10883-MONOMER[c]', 'EG10885-MONOMER[c]', 'EG10890-MONOMER[c]',
										'EG10864-MONOMER[c]', 'EG10881-MONOMER[c]', 'EG10865-MONOMER[c]',
										'EG10868-MONOMER[c]', 'EG10880-MONOMER[c]', 'EG10867-MONOMER[c]',
										'EG10866-MONOMER[c]'],
			's50_proteinComplexes'	:	['CPLX0-3956[c]'],
			's50_23sRRNA'			:	['RRLA-RRNA[c]','RRLB-RRNA[c]','RRLC-RRNA[c]','RRLD-RRNA[c]',
										'RRLE-RRNA[c]','RRLG-RRNA[c]','RRLH-RRNA[c]'],
			's50_5sRRNA'			:	['RRFA-RRNA[c]','RRFB-RRNA[c]','RRFC-RRNA[c]','RRFD-RRNA[c]',
										'RRFE-RRNA[c]','RRFG-RRNA[c]','RRFH-RRNA[c]'],
			's50_fullComplex'		:	['CPLX-50SA[c]'],
		}

		self.__dict__.update(moleculeGroups)