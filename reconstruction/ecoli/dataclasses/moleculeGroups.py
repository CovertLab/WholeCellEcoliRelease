"""
SimulationData molecule groups

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

import reconstruction.ecoli.dataclasses

class moleculeGroups(DataClass.DataClass):
	""" moleculeGroups """

	def __init__(self, simData):
		super(moleculeGroups, self, simData).__init__()
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
		}

		self.__dict__.update(moleculeGroups)