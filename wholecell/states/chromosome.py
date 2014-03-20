
from __future__ import division

import numpy as np
import tables

import wholecell.states.state
from wholecell.containers.chromosome_container import ChromosomeContainer

N_BASES = 5000000 # TODO: from kb
STRAND_MULTIPLICITY = 3 # TODO: estimate somehow? from kb?

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		},
	'DNA polymerase':{
		},
	}


# TODO: design views, queries, requests, etc

class Chromosome(wholecell.states.state.State):

	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'Chromosome',
			'name':'Chromosome',
			'dynamics':[],
			'units':{}
			}

		self.container = None

		super(Chromosome, self).__init__(*args, **kwargs)

	
	def initialize(self, sim, kb):
		super(Chromosome, self).initialize(sim, kb)

		self.container = ChromosomeContainer(N_BASES, STRAND_MULTIPLICITY,
			MOLECULE_ATTRIBUTES)


	def calcInitialConditions(self):
		pass

