"""
SimulationData for Ecoli

Raw data processed into forms convienent for whole-cell modeling

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""
from __future__ import division

from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

class SimulationDataEcoli(object):
	""" SimulationDataEcoli """

	def __init__(self):
		self.raw_data = KnowledgeBaseEcoli()
				