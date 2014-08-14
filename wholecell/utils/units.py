#!/usr/bin/env python

"""
Units

Defines/registers custom units for Pint

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/14/2014
"""

import scipy.constants
from unum.units import *
from unum import Unum

count = Unum.unit('count',mol/(scipy.constants.Avogadro))
nt = Unum.unit('nucleotide', count)
aa = Unum.unit('amino_acid', count)