#!/usr/bin/env python

"""
Logger

Abstract class which defines the interface loggers expose to the simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/28/2013
"""

from __future__ import absolute_import, division, print_function

import abc
import six

class Logger(six.with_metaclass(abc.ABCMeta, object)):
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
