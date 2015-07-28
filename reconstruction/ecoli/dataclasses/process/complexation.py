
from __future__ import division

import numpy as np
from wholecell.utils import units

class Complexation(object):
	def __init__(self, raw_data, sim_data):
		# Build the abstractions needed for complexation

		molecules = []

		subunits = []
		complexes = []

		stoichMatrixI = []
		stoichMatrixJ = []
		stoichMatrixV = []

		stoichMatrixMass = []

		# Remove complexes that are currently not simulated
		FORBIDDEN_MOLECULES = {
			"modified-charged-selC-tRNA", # molecule does not exist
			}

		deleteReactions = []
		for reactionIndex, reaction in enumerate(raw_data.complexationReactions):
			for molecule in reaction["stoichiometry"]:
				if molecule["molecule"] in FORBIDDEN_MOLECULES:
					deleteReactions.append(reactionIndex)
					break

		for reactionIndex in deleteReactions[::-1]:
			del raw_data.complexationReactions[reactionIndex]

		for reactionIndex, reaction in enumerate(raw_data.complexationReactions):
			assert reaction["process"] == "complexation"
			assert reaction["dir"] == 1

			for molecule in reaction["stoichiometry"]:
				if molecule["type"] == "metabolite":
					moleculeName = "{}[{}]".format(
						molecule["molecule"].upper(), # this is stupid # agreed
						molecule["location"]
						)

				else:
					moleculeName = "{}[{}]".format(
						molecule["molecule"],
						molecule["location"]
						)

				if moleculeName not in molecules:
					molecules.append(moleculeName)
					moleculeIndex = len(molecules) - 1

				else:
					moleculeIndex = molecules.index(moleculeName)

				coefficient = molecule["coeff"]

				assert coefficient % 1 == 0


				stoichMatrixI.append(moleculeIndex)
				stoichMatrixJ.append(reactionIndex)
				stoichMatrixV.append(coefficient)

				if coefficient < 0:
					subunits.append(moleculeName)

				else:
					assert molecule["type"] == "proteincomplex"
					complexes.append(moleculeName)

				# Find molecular mass
				molecularMass = sim_data.getter.getMass([moleculeName]).asNumber(units.g / units.mol)[0]
				stoichMatrixMass.append(molecularMass)


		self._stoichMatrixI = np.array(stoichMatrixI)
		self._stoichMatrixJ = np.array(stoichMatrixJ)
		self._stoichMatrixV = np.array(stoichMatrixV)

		self.moleculeNames = molecules
		self.subunitNames = set(subunits)
		self.complexNames = set(complexes)

		# Mass balance matrix
		self._stoichMatrixMass = np.array(stoichMatrixMass)
		self.balanceMatrix = self.stoichMatrix()*self.massMatrix()

		# Find the mass balance of each equation in the balanceMatrix
		massBalanceArray = self.massBalance()

		# The stoichometric matrix should balance out to numerical zero.
		assert np.max([abs(x) for x in massBalanceArray]) < 10e-9

	def stoichMatrix(self):
		shape = (self._stoichMatrixI.max()+1, self._stoichMatrixJ.max()+1)

		out = np.zeros(shape, np.float64)

		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixV

		return out

	def massMatrix(self):
		shape = (self._stoichMatrixI.max()+1, self._stoichMatrixJ.max()+1)

		out = np.zeros(shape, np.float64)

		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixMass

		return out

	def massBalance(self):
		'''
		Sum along the columns of the massBalance matrix to check for reaction
		mass balance
		'''

		reactionSumsArray = []
		
		for index, column in enumerate(self.balanceMatrix.T):
			reactionSumsArray.append(sum(column))

		return reactionSumsArray

	# TODO: redesign this so it doesn't need to create a stoich matrix
	def getMonomers(self, cplxId):
		'''
		Returns subunits for a complex (or any ID passed).
		If the ID passed is already a monomer returns the
		monomer ID again with a stoichiometric coefficient
		of zero.
		'''

		info = self._moleculeRecursiveSearch(cplxId, self.stoichMatrix(), np.array(self.moleculeNames))

		return {'subunitIds' : np.array(info.keys()), 'subunitStoich' : -1 * np.array(info.values())}

	def _findRow(self, product,speciesList):

		for sp in range(0, len(speciesList)):
			if speciesList[sp] == product: return sp
		return -1

	def _findColumn(self, stoichMatrixRow, row):

		for i in range(0,len(stoichMatrixRow)):
			if int(stoichMatrixRow[i]) == 1: return i
		return -1

	def _moleculeRecursiveSearch(self, product, stoichMatrix, speciesList, flag = 0):
		row = self._findRow(product,speciesList)
		if row == -1: return []

		col = self._findColumn(stoichMatrix[row,:], row)
		if col == -1:
			if flag == 0: return []
			else: return {product: -1}

		total = {}
		for i in range(0, len(speciesList)):
			if i == row: continue
			val = stoichMatrix[i][col]
			sp = speciesList[i]

			if val:
				x = self._moleculeRecursiveSearch(sp, stoichMatrix, speciesList, 1)
				for j in x:
					if j in total: total[j] += x[j]*(abs(val))
					else: total[j] = x[j]*(abs(val))
		return total
