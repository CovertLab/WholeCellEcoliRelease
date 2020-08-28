"""
SimulationData for the Complexation process

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 01/23/2015
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from wholecell.utils import units
from six.moves import range, zip


class ComplexationError(Exception):
	pass

class MoleculeNotFoundError(ComplexationError):
	pass

class Complexation(object):
	def __init__(self, raw_data, sim_data):
		# Build the abstractions needed for complexation
		molecules = []  # List of all molecules involved in complexation
		subunits = []  # List of all molecules that participate as subunits
		complexes = []  # List of all molecules that participate as complexes
		stoichMatrixI = []  # Molecule indices
		stoichMatrixJ = []  # Reaction indices
		stoichMatrixV = []  # Stoichiometric coefficients
		stoichMatrixMass = []  # Molecular masses of molecules in stoichMatrixI

		# Remove complexes containing molecules that are currently not simulated
		FORBIDDEN_MOLECULES = {
			"modified-charged-selC-tRNA",  # molecule does not exist
			}

		deleteReactions = []
		for reactionIndex, reaction in enumerate(raw_data.complexationReactions):
			for molecule in reaction["stoichiometry"]:
				if molecule["molecule"] in FORBIDDEN_MOLECULES:
					deleteReactions.append(reactionIndex)
					break

		for reactionIndex in deleteReactions[::-1]:
			del raw_data.complexationReactions[reactionIndex]

		# Build stoichiometric matrix from given complexation reactions
		for reactionIndex, reaction in enumerate(raw_data.complexationReactions):
			assert reaction["process"] == "complexation"
			assert reaction["dir"] == 1

			for molecule in reaction["stoichiometry"]:
				if molecule["type"] == "metabolite":
					moleculeName = "{}[{}]".format(
						molecule["molecule"].upper(),
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
				assert (coefficient % 1) == 0

				stoichMatrixI.append(moleculeIndex)
				stoichMatrixJ.append(reactionIndex)
				stoichMatrixV.append(coefficient)

				# Classify molecule into subunit or complex depending on sign
				# of the stoichiometric coefficient - Note that a molecule can
				# be both a subunit and a complex
				if coefficient < 0:
					subunits.append(moleculeName)
				else:
					assert molecule["type"] == "proteincomplex"
					complexes.append(moleculeName)

				# Find molecular mass of the molecule and add to mass matrix
				molecularMass = sim_data.getter.getMass([moleculeName]).asNumber(units.g / units.mol)[0]
				stoichMatrixMass.append(molecularMass)

		self.rates = np.full(
			(len(raw_data.complexationReactions),),
			raw_data.parameters['complexation_rate'].asNumber())

		self._stoichMatrixI = np.array(stoichMatrixI)
		self._stoichMatrixJ = np.array(stoichMatrixJ)
		self._stoichMatrixV = np.array(stoichMatrixV)
		self._stoichMatrixMass = np.array(stoichMatrixMass)

		self.moleculeNames = molecules
		self.ids_complexes = [self.moleculeNames[i] for i in np.where(np.any(self.stoichMatrix() > 0, axis=1))[0]]
		self.ids_reactions = [stoich_dict['id'] for stoich_dict in raw_data.complexationReactions]

		# Remove duplicate names in subunits and complexes
		self.subunitNames = set(subunits)
		self.complexNames = set(complexes)

		# Create sparse matrix for monomer to complex stoichiometry
		i, j, v, shape = self._buildStoichMatrixMonomers()
		self._stoichMatrixMonomersI = i
		self._stoichMatrixMonomersJ = j
		self._stoichMatrixMonomersV = v
		self._stoichMatrixMonomersShape = shape

		# Mass balance matrix
		# All reaction mass balances should balance out to numerical zero
		balanceMatrix = self.stoichMatrix()*self.massMatrix()
		massBalanceArray = np.sum(balanceMatrix, axis=0)
		assert np.max(np.absolute(massBalanceArray)) < 1e-8  # had to bump this up to 1e-8 because of flagella supercomplex

	def stoichMatrix(self):
		"""
		Builds a stoichiometric matrix based on each given complexation
		reaction. One reaction corresponds to one column in the stoichiometric
		matrix.
		"""
		shape = (self._stoichMatrixI.max() + 1, self._stoichMatrixJ.max() + 1)
		out = np.zeros(shape, np.float64)
		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixV
		return out

	def massMatrix(self):
		"""
		Builds a matrix with the same shape as the stoichiometric matrix, but
		with molecular masses as elements instead of stoichiometric constants
		"""
		shape = (self._stoichMatrixI.max() + 1, self._stoichMatrixJ.max() + 1)
		out = np.zeros(shape, np.float64)
		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixMass
		return out

	def stoichMatrixMonomers(self):
		"""
		Returns the dense stoichiometric matrix for monomers from each complex
		"""
		out = np.zeros(self._stoichMatrixMonomersShape, np.float64)
		out[self._stoichMatrixMonomersI, self._stoichMatrixMonomersJ] = self._stoichMatrixMonomersV
		return out

	# TODO: redesign this so it doesn't need to create a stoich matrix
	def getMonomers(self, cplxId):
		'''
		Returns subunits for a complex (or any ID passed). If the ID passed is
		already a monomer returns the monomer ID again with a stoichiometric
		coefficient of one.
		'''
		info = self._moleculeRecursiveSearch(cplxId, self.stoichMatrix(), self.moleculeNames)
		subunits = []
		subunit_stoich = []
		for subunit, stoich in sorted(info.items()):
			subunits.append(subunit)
			subunit_stoich.append(stoich)
		return {'subunitIds': np.array(subunits), 'subunitStoich': np.array(subunit_stoich)}

	def _buildStoichMatrixMonomers(self):
		"""
		Builds a stoichiometric matrix where each column is a reaction that
		forms a complex directly from its constituent monomers. Since some
		reactions from the raw data are complexation reactions of complexes,
		this is different from the stoichiometric matrix generated by
		stoichMatrix().
		"""
		stoichMatrixMonomersI = []
		stoichMatrixMonomersJ = []
		stoichMatrixMonomersV = []

		for colIdx, id_complex in enumerate(self.ids_complexes):
			D = self.getMonomers(id_complex)

			rowIdx = self.moleculeNames.index(id_complex)
			stoichMatrixMonomersI.append(rowIdx)
			stoichMatrixMonomersJ.append(colIdx)
			stoichMatrixMonomersV.append(1.)

			for subunitId, subunitStoich in zip(D["subunitIds"], D["subunitStoich"]):
				rowIdx = self.moleculeNames.index(subunitId)
				stoichMatrixMonomersI.append(rowIdx)
				stoichMatrixMonomersJ.append(colIdx)
				stoichMatrixMonomersV.append(-1. * subunitStoich)

		stoichMatrixMonomersI = np.array(stoichMatrixMonomersI)
		stoichMatrixMonomersJ = np.array(stoichMatrixMonomersJ)
		stoichMatrixMonomersV = np.array(stoichMatrixMonomersV)

		shape = (stoichMatrixMonomersI.max() + 1, stoichMatrixMonomersJ.max() + 1)

		return (stoichMatrixMonomersI, stoichMatrixMonomersJ, stoichMatrixMonomersV, shape)

	def _findRow(self, product, speciesList):
		try:
			row = speciesList.index(product)
		except ValueError as e:
			raise MoleculeNotFoundError(
				"Could not find %s in the list of molecules." % (product,), e)
		return row

	def _findColumn(self, stoichMatrixRow):
		for i in range(0, len(stoichMatrixRow)):
			if int(stoichMatrixRow[i]) == 1:
				return i
		return -1  # Flag for monomer

	def _moleculeRecursiveSearch(self, product, stoichMatrix, speciesList):
		row = self._findRow(product, speciesList)
		col = self._findColumn(stoichMatrix[row, :])
		if col == -1:
			return {product: 1.0}

		total = {}
		for i in range(0, len(speciesList)):
			if i == row:
				continue
			val = stoichMatrix[i][col]
			sp = speciesList[i]

			if val != 0:
				x = self._moleculeRecursiveSearch(sp, stoichMatrix, speciesList)
				for j in x:
					if j in total:
						total[j] += x[j]*(np.absolute(val))
					else:
						total[j] = x[j]*(np.absolute(val))
		return total
