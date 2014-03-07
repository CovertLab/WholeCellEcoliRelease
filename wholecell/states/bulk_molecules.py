#!/usr/bin/env python

"""
BulkMolecules.py

State which represents for a class of molecules the bulk copy numbers.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/04/2013
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

# TODO: break this into multiple files, it's becoming unbearably long

import re

import numpy as np
import tables

import wholecell.states.state
import wholecell.states.partition
import wholecell.utils.bulk_objects_container


class BulkMolecules(wholecell.states.state.State):
	pass


class BulkMoleculesPartition(wholecell.states.partition.Partition):
	pass
