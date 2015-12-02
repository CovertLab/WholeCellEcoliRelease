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

	def _initialize(self, validation_data_raw):

		self._loadProteinDatasets(validation_data_raw)


	def _loadProteinDatasets(self, validation_data_raw):

		# Load taniguichi Xie Science 2010 dataset
		taniguichi_dataset = validation_data_raw.taniguichi2010_table_6
		self.taniguichi2010counts = np.zeros(len(taniguichi_dataset), dtype=[('gene_symbol', '|S10'), ('b_number', '|S10'), ('counts_ave', np.float32), ('gamma_shape_parameter', np.float32), ('gamma_scale_parameter', np.float32)])
		for idx, row in enumerate(taniguichi_dataset):
			self.taniguichi2010counts[idx]["gene_symbol"] = row["Gene Name"]
			self.taniguichi2010counts[idx]["b_number"] = row["B Number"]
			self.taniguichi2010counts[idx]["counts_ave"] = row["Mean_Protein"]
			self.taniguichi2010counts[idx]["gamma_shape_parameter"] = row["A_Protein"]
			self.taniguichi2010counts[idx]["gamma_scale_parameter"] = row["B_Protein"]

		# Load Houser Wilke PLoSCB 2015 dataset
		houser_dataset = validation_data_raw.houser2015_javier_table
		self.houser2015counts = np.zeros(len(houser_dataset), dtype=[
			('gene_symbol', '|S10'),
			('counts_ave_exponential', np.float32),
			('sample16_t3', np.float32),
			('sample17_t4', np.float32),
			('sample18_t5', np.float32),
			('sample19_t6', np.float32),
			('sample20_t8', np.float32),
			('sample21_t24', np.float32),
			('sample22_t48', np.float32),
			('sample23_t168', np.float32),
			('sample24_336', np.float32),
			('sample25_t3', np.float32),
			('sample26_t4', np.float32),
			('sample27_t5', np.float32),
			('sample28_t6', np.float32),
			('sample29_t8', np.float32),
			('sample30_t24', np.float32),
			('sample31_t48', np.float32),
			('sample32_t168', np.float32),
			('sample34_336', np.float32),
			('sample97_t3', np.float32),
			('sample98_t4', np.float32),
			('sample99_t5', np.float32),
			('sample100_t6', np.float32),
			('sample101_t8', np.float32),
			('sample102_t24', np.float32),
			('sample103_t48', np.float32),
			('sample104_t168', np.float32),
			('sample105_336', np.float32)
			])

		for idx, row in enumerate(houser_dataset):
			self.houser2015counts[idx]["gene_symbol"] = row["gene_symbol"]
			# Hour 6 is the most clearly exponential time point, this is in columns labeled as t6
			exponential_average = np.mean([float(row["sample19_t6"]), float(row["sample28_t6"]), float(row["sample100_t6"])])
			# Only include proteins observed at least once in this time step
			if exponential_average:
				self.houser2015counts[idx]["counts_ave_exponential"] = exponential_average
			# Otherwise, remove this row from the dataset
			else:
				np.delete(self.houser2015counts, idx, 0)

			# Load the rest of the data as-is
			for fieldName in row:
				if fieldName == "gene_symbol" or fieldName == "":
					continue
				self.houser2015counts[idx][fieldName] = row[fieldName]