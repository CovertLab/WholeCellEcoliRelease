"""
SimulationData for the Complexation process
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

		# Get IDs of reactions that should be removed
		removed_reaction_ids = {
			rxn['id'] for rxn in raw_data.complexation_reactions_removed}

		self.ids_reactions = []
		reaction_index = 0

		# Build stoichiometric matrix from given complexation reactions
		for reaction in raw_data.complexation_reactions:
			# Skip removed reactions
			if reaction['id'] in removed_reaction_ids:
				continue

			self.ids_reactions.append(reaction['id'])

			for mol_id, coeff in reaction["stoichiometry"].items():
				mol_id_with_compartment = "{}[{}]".format(
					mol_id,
					sim_data.getter.get_compartment(mol_id)[0]
					)

				if mol_id_with_compartment not in molecules:
					molecules.append(mol_id_with_compartment)
					molecule_index = len(molecules) - 1
				else:
					molecule_index = molecules.index(mol_id_with_compartment)

				# Assume coefficents given as null are -1
				if coeff is None:
					coeff = -1

				assert (coeff % 1) == 0

				stoichMatrixI.append(molecule_index)
				stoichMatrixJ.append(reaction_index)
				stoichMatrixV.append(coeff)

				# Classify molecule into subunit or complex depending on sign
				# of the stoichiometric coefficient - Note that a molecule can
				# be both a subunit and a complex
				if coeff < 0:
					subunits.append(mol_id_with_compartment)
				else:
					complexes.append(mol_id_with_compartment)

				# Find molecular mass of the molecule and add to mass matrix
				molecularMass = sim_data.getter.get_mass(mol_id_with_compartment).asNumber(units.g / units.mol)
				stoichMatrixMass.append(molecularMass)

			reaction_index += 1

		self.rates = np.full((reaction_index, ),
			sim_data.constants.complexation_rate.asNumber(1/units.s))

		self._stoich_matrix_I = np.array(stoichMatrixI)
		self._stoich_matrix_J = np.array(stoichMatrixJ)
		self._stoich_matrix_V = np.array(stoichMatrixV)
		self._stoich_matrix_mass = np.array(stoichMatrixMass)

		self.molecule_names = molecules
		self.ids_complexes = [self.molecule_names[i] for i in np.where(np.any(self.stoich_matrix() > 0, axis=1))[0]]

		# Remove duplicate names in subunits and complexes
		self.subunit_names = set(subunits)
		self.complex_names = set(complexes)

		# Create sparse matrix for monomer to complex stoichiometry
		i, j, v, shape = self._buildStoichMatrixMonomers()
		self._stoichMatrixMonomersI = i
		self._stoichMatrixMonomersJ = j
		self._stoichMatrixMonomersV = v
		self._stoichMatrixMonomersShape = shape

		# Mass balance matrix
		# All reaction mass balances should balance out to numerical zero
		balanceMatrix = self.stoich_matrix() * self.mass_matrix()
		massBalanceArray = np.sum(balanceMatrix, axis=0)
		assert np.max(np.absolute(massBalanceArray)) < 1e-8  # had to bump this up to 1e-8 because of flagella supercomplex

	def stoich_matrix(self):
		"""
		Builds a stoichiometric matrix based on each given complexation
		reaction. One reaction corresponds to one column in the stoichiometric
		matrix.
		"""
		shape = (self._stoich_matrix_I.max() + 1, self._stoich_matrix_J.max() + 1)
		out = np.zeros(shape, np.float64)
		out[self._stoich_matrix_I, self._stoich_matrix_J] = self._stoich_matrix_V
		return out

	def mass_matrix(self):
		"""
		Builds a matrix with the same shape as the stoichiometric matrix, but
		with molecular masses as elements instead of stoichiometric constants
		"""
		shape = (self._stoich_matrix_I.max() + 1, self._stoich_matrix_J.max() + 1)
		out = np.zeros(shape, np.float64)
		out[self._stoich_matrix_I, self._stoich_matrix_J] = self._stoich_matrix_mass
		return out

	def stoich_matrix_monomers(self):
		"""
		Returns the dense stoichiometric matrix for monomers from each complex
		"""
		out = np.zeros(self._stoichMatrixMonomersShape, np.float64)
		out[self._stoichMatrixMonomersI, self._stoichMatrixMonomersJ] = self._stoichMatrixMonomersV
		return out

	# TODO: redesign this so it doesn't need to create a stoich matrix
	def get_monomers(self, cplxId):
		'''
		Returns subunits for a complex (or any ID passed). If the ID passed is
		already a monomer returns the monomer ID again with a stoichiometric
		coefficient of one.
		'''
		info = self._moleculeRecursiveSearch(cplxId, self.stoich_matrix(), self.molecule_names)
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
			D = self.get_monomers(id_complex)

			rowIdx = self.molecule_names.index(id_complex)
			stoichMatrixMonomersI.append(rowIdx)
			stoichMatrixMonomersJ.append(colIdx)
			stoichMatrixMonomersV.append(1.)

			for subunitId, subunitStoich in zip(D["subunitIds"], D["subunitStoich"]):
				rowIdx = self.molecule_names.index(subunitId)
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
