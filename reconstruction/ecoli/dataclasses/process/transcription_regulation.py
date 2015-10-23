
from __future__ import division

import numpy as np
import os
import theano.tensor as T
import theano
import cPickle
import wholecell
from wholecell.utils import units
from . import metabolism

class TranscriptionRegulation(object):
	def __init__(self, raw_data, sim_data):
		return