#!/usr/bin/env python

"""
Fitter

Adjust simulation parameters

TODO:
- document the math
- compute and use activation rates for RNA poly, ribosomes
- fit metabolism enzyme expression
- replace fake metabolite targets with measured metabolite targets

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/11/2013
"""

from __future__ import division

import numpy as np
import os
import collections

from wholecell.utils import units

from reconstruction.ecoli.fitkb1 import fitKb_1
from reconstruction.ecoli.fitkb2 import fitKb_2

def fitAtLevel(fitLevel, kb, simOutDir):
	# TODO: Obviously make this more sophisticated
	if fitLevel == 1:
		fitKb_1(kb)

	if fitLevel == 2:
		fitKb_2(kb, simOutDir)