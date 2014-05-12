#!/usr/bin/env python

"""
Test replication.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 4/6/2014
"""
import wholecell.utils.polymerize as p

import numpy as np

import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest

class Test_polymerize(unittest.TestCase):
	
	@noseAttrib.attr('polymerizeTest')
	@noseAttrib.attr('smalltest')
	def test_polymerize_functionCall(self):
		sequences = np.matrix([['A','T','A','A'],['G','T','A','C'],['C','C','T','T']])
		baseAmounts = np.array([9,9,9,9])
		bases = np.array(['A','C','G','T'])
		basePadValue = ' '
		energy = 0
		energyCostPerBase = 0
		
		# Good calls test
		p.polymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
	