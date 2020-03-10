"""
SimulationData process associated data

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/12/2015
"""

from __future__ import absolute_import, division, print_function

from .chromosome_structure import ChromosomeStructure
from .complexation import Complexation
from .equilibrium import Equilibrium
from .metabolism import Metabolism
from .replication import Replication
from .rna_decay import RnaDecay
from .transcription import Transcription
from .transcription_regulation import TranscriptionRegulation
from .translation import Translation
from .two_component_system import TwoComponentSystem


class Process(object):
	""" Process """

	def __init__(self, raw_data, sim_data):

		self.chromosome_structure = ChromosomeStructure(raw_data, sim_data)
		self.complexation = Complexation(raw_data, sim_data)
		self.equilibrium = Equilibrium(raw_data, sim_data)
		self.metabolism = Metabolism(raw_data, sim_data)
		self.replication = Replication(raw_data, sim_data)
		self.rna_decay = RnaDecay(raw_data, sim_data)
		self.transcription = Transcription(raw_data, sim_data)
		self.transcription_regulation = TranscriptionRegulation(raw_data, sim_data)
		self.translation = Translation(raw_data, sim_data)
		self.two_component_system = TwoComponentSystem(raw_data, sim_data)
