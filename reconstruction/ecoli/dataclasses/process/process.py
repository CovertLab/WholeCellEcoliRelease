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
from reconstruction.ecoli.dataclasses.process.complexation import Complexation
from reconstruction.ecoli.dataclasses.process.rna_decay import RnaDecay
from reconstruction.ecoli.dataclasses.process.equilibrium import Equilibrium
from reconstruction.ecoli.dataclasses.process.transcription_regulation import TranscriptionRegulation
from reconstruction.ecoli.dataclasses.process.two_component_system import TwoComponentSystem


import re
import numpy as np

class Process(object):
	""" Process """

	def __init__(self, raw_data, sim_data):

		self.replication = Replication(raw_data, sim_data)
		self.metabolism = Metabolism(raw_data, sim_data)
		self.transcription = Transcription(raw_data, sim_data)
		self.translation = Translation(raw_data, sim_data)
		self.complexation = Complexation(raw_data, sim_data)
		self.rna_decay = RnaDecay(raw_data, sim_data)
		self.equilibrium = Equilibrium(raw_data, sim_data)
		self.transcription_regulation = TranscriptionRegulation(raw_data, sim_data)
		self.two_component_system = TwoComponentSystem(raw_data, sim_data)
