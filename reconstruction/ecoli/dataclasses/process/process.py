"""
SimulationData process associated data

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import division

from reconstruction.ecoli.dataclasses.process.replication import Replication
from reconstruction.ecoli.dataclasses.process.metabolism import Metabolism
from reconstruction.ecoli.dataclasses.process.transcription import Transcription
from reconstruction.ecoli.dataclasses.process.translation import Translation



import re
import numpy as np

class Process(object):
	""" Process """

	def __init__(self, raw_data, sim_data):

		self.replication = Replication(raw_data, sim_data)
		self.metabolism = Metabolism(raw_data, sim_data)
		self.transcription = Transcription(raw_data, sim_data)
		self.translation = Translation(raw_data, sim_data)