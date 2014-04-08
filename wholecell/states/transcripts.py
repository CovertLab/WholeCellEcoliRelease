
from __future__ import division

# import numpy as np
# import tables

import wholecell.states.state
from wholecell.containers.transcripts_container import TranscriptsContainer

ARRAY_LENGTH = 1000000 # TODO: estimate somehow?

MOLECULE_ATTRIBUTES = {
	}


# TODO: design views, queries, requests, etc

class Transcripts(wholecell.states.state.State):

	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'Transcripts',
			'name':'Transcripts',
			'dynamics':[],
			'units':{}
			}

		self.container = None

		super(Transcripts, self).__init__(*args, **kwargs)

	
	def initialize(self, sim, kb):
		super(Transcripts, self).initialize(sim, kb)

		self.container = TranscriptsContainer(ARRAY_LENGTH,
			MOLECULE_ATTRIBUTES, self.randStream)


	def calcInitialConditions(self):
		pass


	def pytablesCreate(self, h5file, expectedRows):
		self.container.pytablesCreate(h5file)


	def pytablesAppend(self, h5file):
		self.container.pytablesAppend(h5file, self.time.value)


	def pytablesLoad(self, h5file, timePoint):
		self.container.pytablesLoad(h5file, timePoint)

