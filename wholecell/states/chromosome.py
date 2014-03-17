
from __future__ import division

import numpy as np
import tables

import wholecell.states.state
import wholecell.utils.chromosome_bound_molecules_container

N_BASES = 5000000 # TODO: from kb
STRAND_MULTIPLICITY = 3
MOLECULE_WIDTH = 50 # TODO: from kb

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

		# self.time = None

		self.container = None

		super(Chromosome, self).__init__(*args, **kwargs)

	
	def initialize(self, sim, kb):
		super(Chromosome, self).initialize(sim, kb)

		self.container = wholecell.utils.chromosome_bound_molecules_container.ChromosomeBoundMoleculeContainer(
			N_BASES, STRAND_MULTIPLICITY, MOLECULE_ATTRIBUTES)


	def calcInitialConditions(self):
		# Add replication fork

		forkStartPos = 20000
		forkStopPos = 30000 - 1

		(forkStart, forkStop) = self.container.divideRegion('R', forkStartPos, forkStopPos)

		polyStart = self.container.moleculeNew('DNA polymerase')
		polyStop = self.container.moleculeNew('DNA polymerase')

		forwardExtent = 50
		reverseExtent = 50

		nSteps = 3600

		ntRoot = np.zeros(nSteps, np.int64)
		ntChildA = np.zeros(nSteps, np.int64)
		ntChildB = np.zeros(nSteps, np.int64)

		try:
			for i in xrange(nSteps):
				self.container.moleculeLocationIsUnbound(polyStart)
				self.container.moleculeLocationIsUnbound(polyStop)

				self.container.forkExtend(forkStart, 100)
				self.container.forkExtend(forkStop, 100)

				self.container.moleculeLocationIsFork(polyStart, forkStart, forwardExtent, reverseExtent)
				self.container.moleculeLocationIsFork(polyStop, forkStop, forwardExtent, reverseExtent)

				ntRoot[i] = (self.container._array[0, :] != self.container._inactive).sum()
				ntChildA[i] = (self.container._array[1, :] != self.container._inactive).sum()
				ntChildB[i] = (self.container._array[2, :] != self.container._inactive).sum()


		except Exception as e:
			print e

		import matplotlib.pyplot as plt

		plt.plot(ntRoot[:i])
		plt.plot(ntChildA[:i])
		plt.plot(ntChildB[:i])

		plt.show()

		import ipdb; ipdb.set_trace()

		# assert self.container.moleculeLocation(molecule) is None

		# self.container.moleculeLocationIs(
		# 	molecule,
		# 	'R', 100, '+',
		# 	50, 50
		# 	)

		# self.container.moleculeLocationIs(
		# 	molecule,
		# 	'R', 200, '+',
		# 	50, 50
		# 	)

		# assert self.container.moleculeLocation(molecule) == ('R', 100, '+', 50, 50)

		# assert self.container.moleculesBound() == {molecule}
		# assert self.container.moleculesBound('DNA polymerase') == {molecule}
		# assert self.container.moleculesBound('RNA polymerase') == set()

		# assert self.container.moleculesBound(
		# 	strand = 'R',
		# 	position = 100,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	strand = 'R',
		# 	position = 140,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	strand = 'R',
		# 	position = 160,
		# 	) == set()

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 100,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	moleculeName = 'RNA polymerase',
		# 	strand = 'R',
		# 	position = 100,
		# 	) == set()

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 100,
		# 	direction = '+',
		# 	extentForward = 10,
		# 	extentReverse = 10,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 100,
		# 	direction = '-',
		# 	extentForward = 10,
		# 	extentReverse = 10,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 160,
		# 	direction = '-',
		# 	extentForward = 20,
		# 	extentReverse = 0,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 160,
		# 	direction = '+',
		# 	extentForward = 20,
		# 	extentReverse = 0,
		# 	) == set()

		# import ipdb; ipdb.set_trace()

