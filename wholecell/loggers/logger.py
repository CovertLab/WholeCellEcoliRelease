"""
Logger

Abstract class which defines the interface loggers expose to the simulation
"""

import abc


class Logger(metaclass=abc.ABCMeta):
	""" Logger """

	@abc.abstractmethod
	def initialize(self, sim):
		""" initialize -- called at beginning of simulation """
		return

	@abc.abstractmethod
	def append(self, sim):
		""" append -- called at each iteration of simulation """
		return

	@abc.abstractmethod
	def finalize(self, sim):
		""" finalize -- called at end of simulation """
		return
