"""
SimulationData for Ecoli

Raw data processed into forms convienent for whole-cell modeling

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""
from __future__ import division

import numpy as np
import collections

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

# Data classes
from reconstruction.ecoli.dataclasses.getterFunctions import getterFunctions
from reconstruction.ecoli.dataclasses.moleculeGroups import moleculeGroups
from reconstruction.ecoli.dataclasses.state.state import State
from reconstruction.ecoli.dataclasses.process.process import Process

class SimulationDataEcoli(object):
	""" SimulationDataEcoli """

	def __init__(self):
		# Raw data
		raw_data = KnowledgeBaseEcoli()
		self._addHardCodedAttributes()

		# Data classes
		self.getter = getterFunctions(raw_data, self)
		self.moleculeGroups = moleculeGroups(raw_data, self)
		self.process = Process(raw_data, self)
		self.state = State(raw_data, self)

	def _addHardCodedAttributes(self):
		self.molecular_weight_keys = [
			'23srRNA',
			'16srRNA',
			'5srRNA',
			'tRNA',
			'mRNA',
			'miscRNA',
			'protein',
			'metabolite',
			'water',
			'DNA',
			'RNA' # nonspecific RNA
			]

		self.molecular_weight_order = collections.OrderedDict([
			(key, index) for index, key in enumerate(self.molecular_weight_keys)
			])

		self.amino_acid_1_to_3_ordered = collections.OrderedDict((
			("A", "ALA-L[c]"), ("R", "ARG-L[c]"), ("N", "ASN-L[c]"), ("D", "ASP-L[c]"),
			("C", "CYS-L[c]"), ("E", "GLU-L[c]"), ("Q", "GLN-L[c]"), ("G", "GLY[c]"),
			("H", "HIS-L[c]"), ("I", "ILE-L[c]"), ("L", "LEU-L[c]"), ("K", "LYS-L[c]"),
			("M", "MET-L[c]"), ("F", "PHE-L[c]"), ("P", "PRO-L[c]"), ("S", "SER-L[c]"),
			("T", "THR-L[c]"), ("W", "TRP-L[c]"), ("Y", "TYR-L[c]"), ("U", "SEC-L[c]"),
			("V", "VAL-L[c]")
			))