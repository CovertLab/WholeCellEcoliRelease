#!/usr/bin/env python

"""
UnitRegistration

Defines/registers custom units for Pint

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/9/2014
"""

from pint import UnitRegistry
UREG = UnitRegistry()
Q_ = UREG.Quantity
UREG.define('nucleotide = []')
UREG.define('amino_acid = []')
UREG.define('DCW- = 1.')