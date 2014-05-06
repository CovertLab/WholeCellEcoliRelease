#!/usr/bin/env python

"""
Polymerize

Wrapper around C extension to polymerize

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/06/14
"""

import numpy as np
import wholecell.utils._polymerize

def polymerize(elngRate, deficitAAs, requiredAAs, aaCounts, updatedAAs, aasUsed, seed):

	if type(elngRate) != float:
		raise TypeError, "elngRate must be of type float"

	if type(deficitAAs) != np.ndarray or deficitAAs.dtype != np.int64:
		raise TypeError, "deficitAAs must be a numpy array with dtype np.int64"

	if type(requiredAAs) != np.ndarray or requiredAAs.dtype != np.int64:
		raise TypeError, "requiredAAs must be a numpy array with dtype np.int64"

	if type(aaCounts) != np.ndarray or aaCounts.dtype != np.int64:
		raise TypeError, "aaCounts must be a numpy array with dtype np.int64"

	if type(updatedAAs) != np.ndarray or updatedAAs.dtype != np.int64:
		raise TypeError, "updatedAAs must be a numpy array with dtype np.int64"

	if type(aasUsed) != np.ndarray or aasUsed.dtype != np.int64:
		raise TypeError, "aasUsed must be a numpy array with dtype np.int64"

	if type(seed) != np.uint32:
		raise TypeError, "seed must be of type np.uint32"

	wholecell.utils._polymerize.polymerize(
		elngRate, deficitAAs, requiredAAs, aaCounts,
		updatedAAs, aasUsed, seed
		)