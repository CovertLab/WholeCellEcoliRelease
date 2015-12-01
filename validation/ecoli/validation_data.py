"""
ValidationData for Ecoli

Raw data processed into forms convienent for validation and analysis

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 11/30/2015
"""
from __future__ import division

import numpy as np
import collections
from unum import Unum

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

# Data classes
from reconstruction.ecoli.dataclasses.getterFunctions import getterFunctions
from reconstruction.ecoli.dataclasses.moleculeGroups import moleculeGroups
from reconstruction.ecoli.dataclasses.constants import Constants
from reconstruction.ecoli.dataclasses.state.state import State
from reconstruction.ecoli.dataclasses.process.process import Process
from reconstruction.ecoli.dataclasses.growthRateDependentParameters import Mass, GrowthRateParameters
from reconstruction.ecoli.dataclasses.relation import Relation

class ValidationDataEcoli(object):
	""" ValidationDataEcoli """

	def __init__(self):
		pass

	def initialize(self, validation_data_raw):

		self._loadProteinDatasets(validation_data_raw)


	def _loadProteinDatasets(self, validation_data_raw):

		# Load taniguichi Xie Science 2010 dataset
		taniguichi_dataset = validation_data_raw.taniguichi2010_table_6
		self.taniguichi2010counts = np.zeros(len(taniguichi_dataset), dtype=[('gene_symbol', '|S10'), ('b_number', '|S10'), ('counts_ave', np.float32), ('gamma_shape_parameter', np.float32), ('gamma_scale_parameter', np.float32)])
		for idx, row in enumerate(validation_data_raw.taniguichi2010_table_6):
			self.taniguichi2010counts[idx]["gene_symbol"] = row["Gene Name"]
			self.taniguichi2010counts[idx]["b_number"] = row["B Number"]
			self.taniguichi2010counts[idx]["counts_ave"] = row["Mean_Protein"]
			self.taniguichi2010counts[idx]["gamma_shape_parameter"] = row["A_Protein"]
			self.taniguichi2010counts[idx]["gamma_scale_parameter"] = row["B_Protein"]
