"""
Equilibrium.

TODOs:
fluxesAndMoleculesToSS()
	Consider relocating (since it's useful for both the parca and simulation)
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from scipy import integrate
import sympy as sp

from wholecell.utils import build_ode, data, units
from wholecell.utils.random import stochasticRound
from six.moves import range, zip


class EquilibriumError(Exception):
	pass

class MoleculeNotFoundError(EquilibriumError):
	pass

class Equilibrium(object):
	def __init__(self, raw_data, sim_data):
		# Build the abstractions needed for complexation
		molecules = []

		ratesFwd = []
		ratesRev = []
		rxnIds = []

		stoichMatrixI = []
		stoichMatrixJ = []
		stoichMatrixV = []

		stoichMatrixMass = []
		self.metabolite_set = set()
		self.complex_name_to_rxn_idx = {}

		# Make sure reactions are not duplicated in complexationReactions and equilibriumReactions
		equilibriumReactionIds = set([x["id"] for x in raw_data.equilibrium_reactions])
		complexationReactionIds = set([x["id"] for x in raw_data.complexation_reactions])

		if equilibriumReactionIds.intersection(complexationReactionIds) != set():
			raise Exception(
				"The following reaction ids are specified in equilibriumReactions and complexationReactions: %s" % (
					equilibriumReactionIds.intersection(
						complexationReactionIds)))

		# Remove complexes that are currently not simulated
		FORBIDDEN_MOLECULES = {
			"modified-charged-selC-tRNA", # molecule does not exist
			}

		# Remove reactions that we know won't occur (e.g., don't do computations on metabolites that have zero counts)
		MOLECULES_THAT_WILL_EXIST_IN_SIMULATION = [m["Metabolite"] for m in raw_data.metabolite_concentrations] + ["LEU", "S-ADENOSYLMETHIONINE", "ARABINOSE", "4FE-4S"] + [l["molecules"]["LIGAND"] for l in raw_data.two_component_systems]

		deleteReactions = []
		for reactionIndex, reaction in enumerate(raw_data.equilibrium_reactions):
			for molecule in reaction["stoichiometry"]:
				if molecule["molecule"] in FORBIDDEN_MOLECULES or (molecule["type"] == "metabolite" and molecule["molecule"] not in MOLECULES_THAT_WILL_EXIST_IN_SIMULATION):
					deleteReactions.append(reactionIndex)
					break

		for reactionIndex in deleteReactions[::-1]:
			del raw_data.equilibrium_reactions[reactionIndex]

		# Build stoichiometry matrix
		for reactionIndex, reaction in enumerate(raw_data.equilibrium_reactions):
			assert reaction["process"] == "equilibrium"
			assert reaction["dir"] == 1

			ratesFwd.append(reaction["forward rate"])
			ratesRev.append(reaction["reverse rate"])
			rxnIds.append(reaction["id"])

			for molecule in reaction["stoichiometry"]:
				if molecule["type"] == "metabolite":
					moleculeName = "{}[{}]".format(
						molecule["molecule"].upper(),
						molecule["location"]
						)
					self.metabolite_set.add(moleculeName)

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

				# Store indices for the row and column, and molecule coefficient for building the stoichiometry matrix
				stoichMatrixI.append(moleculeIndex)
				stoichMatrixJ.append(reactionIndex)
				stoichMatrixV.append(coefficient)

				if coefficient > 0:
					assert molecule["type"] == "proteincomplex"
					self.complex_name_to_rxn_idx[moleculeName] = reactionIndex

				# Find molecular mass
				molecularMass = sim_data.getter.get_mass([moleculeName]).asNumber(units.g / units.mol)[0]
				stoichMatrixMass.append(molecularMass)

		# TODO(jerry): Move the rest to a subroutine for __init__ and __setstate__?
		self._stoichMatrixI = np.array(stoichMatrixI)
		self._stoichMatrixJ = np.array(stoichMatrixJ)
		self._stoichMatrixV = np.array(stoichMatrixV)

		self.molecule_names = molecules
		self.ids_complexes = [self.molecule_names[i] for i in np.where(np.any(self.stoich_matrix() > 0, axis=1))[0]]
		self.rxn_ids = rxnIds
		self.rates_fwd = np.array(ratesFwd)
		self.rates_rev = np.array(ratesRev)

		# Mass balance matrix
		self._stoichMatrixMass = np.array(stoichMatrixMass)
		self.balance_matrix = self.stoich_matrix() * self.mass_matrix()

		# Find the mass balance of each equation in the balanceMatrix
		massBalanceArray = self.mass_balance()

		# The stoichometric matrix should balance out to numerical zero.
		assert np.max(np.absolute(massBalanceArray)) < 1e-9

		# Build matrices
		self._populateDerivativeAndJacobian()
		self._stoichMatrix = self.stoich_matrix()

	def __getstate__(self):
		"""Return the state to pickle, omitting derived attributes that
		__setstate__() will recompute, esp. those like the rates for ODEs
		that don't pickle.
		"""
		return data.dissoc_strict(self.__dict__, (
			'_stoichMatrix',
			'Rp', 'Pp', 'mets_to_rxn_fluxes',
			'symbolic_rates', 'symbolic_rates_jacobian',
			'_rates', '_rates_jacobian'))

	def __setstate__(self, state):
		"""Restore instance attributes, recomputing some of them."""
		self.__dict__.update(state)
		self._stoichMatrix = self.stoich_matrix()
		self._populateDerivativeAndJacobian()

	def stoich_matrix(self):
		'''
		Builds stoichiometry matrix
		Rows: molecules
		Columns: reactions
		Values: reaction stoichiometry
		'''
		shape = (self._stoichMatrixI.max()+1, self._stoichMatrixJ.max()+1)
		out = np.zeros(shape, np.float64)
		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixV
		return out

	def mass_matrix(self):
		'''
		Builds stoichiometry mass matrix
		Rows: molecules
		Columns: reactions
		Values: molecular mass
		'''
		shape = (self._stoichMatrixI.max()+1, self._stoichMatrixJ.max()+1)
		out = np.zeros(shape, np.float64)
		out[self._stoichMatrixI, self._stoichMatrixJ] = self._stoichMatrixMass
		return out

	def mass_balance(self):
		'''
		Sum along the columns of the massBalance matrix to check for reaction
		mass balance
		'''
		return np.sum(self.balance_matrix, axis=0)

	def stoich_matrix_monomers(self):
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

		out = np.zeros(shape, np.float64)
		out[stoichMatrixMonomersI, stoichMatrixMonomersJ] = stoichMatrixMonomersV
		return out

	def _populateDerivativeAndJacobian(self):
		'''Compile callable functions for computing the derivative and the Jacobian.'''
		self._makeMatrices()
		self._make_rates()

		self._rates = build_ode.rates(self.symbolic_rates)
		self._rates_jacobian = build_ode.rates_jacobian(
			self.symbolic_rates_jacobian)

	def _makeMatrices(self):
		'''
		Creates matrix that maps metabolites to the flux through their reactions.
		'''
		EPS = 1e-9

		S = self.stoich_matrix()
		Rp = -1. * (S < -1 * EPS) * S
		Pp =  1. * (S >  1 * EPS) * S
		self.Rp = Rp
		self.Pp = Pp

		mets_to_rxn_fluxes = self.stoich_matrix().copy()

		mets_to_rxn_fluxes[(np.abs(mets_to_rxn_fluxes) > EPS).sum(axis = 1) > 1, : ] = 0
		for colIdx in range(mets_to_rxn_fluxes.shape[1]):
			try:
				firstNonZeroIdx = np.where(np.abs(mets_to_rxn_fluxes[:, colIdx]) > EPS)[0][0]
			except IndexError:
				raise Exception(
					"Column %d of S matrix not linearly independent!" % colIdx)
			mets_to_rxn_fluxes[:firstNonZeroIdx, colIdx] = 0
			mets_to_rxn_fluxes[(firstNonZeroIdx + 1):, colIdx] = 0

		self.mets_to_rxn_fluxes = mets_to_rxn_fluxes.T

	def _make_rates(self):
		'''
		Creates symbolic representation of the rates for ordinary differential
		equations and the Jacobian. Used during simulations.
		'''
		S = self.stoich_matrix()

		yStrings = ["y[%d]" % x for x in range(S.shape[0])]
		ratesFwdStrings = ["kf[%d]" % x for x in range(S.shape[0])]
		ratesRevStrings = ["kr[%d]" % x for x in range(S.shape[0])]
		y = sp.symbols(yStrings)
		ratesFwd = sp.symbols(ratesFwdStrings)
		ratesRev = sp.symbols(ratesRevStrings)

		rates = []

		for colIdx in range(S.shape[1]):
			negIdxs = np.where(S[:, colIdx] < 0)[0]
			posIdxs = np.where(S[:, colIdx] > 0)[0]

			reactantFlux = ratesFwd[colIdx]
			for negIdx in negIdxs:
				stoich = -S[negIdx, colIdx]
				if stoich == 1:
					reactantFlux *= y[negIdx]
				else:
					reactantFlux *= y[negIdx]**stoich

			productFlux = ratesRev[colIdx]
			for posIdx in posIdxs:
				stoich = S[posIdx, colIdx]
				if stoich == 1:
					productFlux *= y[posIdx]
				else:
					productFlux *= y[posIdx]**stoich

			rates.append(reactantFlux - productFlux)

		dy = sp.Matrix(rates)
		J = dy.jacobian(y)

		self.symbolic_rates = dy
		self.symbolic_rates_jacobian = J

	def derivatives(self, t, y):
		return self._stoichMatrix.dot(self._rates[0](t, y, self.rates_fwd, self.rates_rev))

	def derivatives_jacobian(self, t, y):
		return self._stoichMatrix.dot(self._rates_jacobian[0](t, y, self.rates_fwd, self.rates_rev))

	def derivatives_jit(self, t, y):
		return self._stoichMatrix.dot(self._rates[1](t, y, self.rates_fwd, self.rates_rev))

	def derivatives_jacobian_jit(self, t, y):
		return self._stoichMatrix.dot(self._rates_jacobian[1](t, y, self.rates_fwd, self.rates_rev))

	def fluxes_and_molecules_to_SS(self, moleculeCounts, cellVolume, nAvogadro,
			random_state, time_limit=1e20, max_iter=100, jit=True):
		y_init = moleculeCounts / (cellVolume * nAvogadro)

		# In this version of SciPy, solve_ivp does not support args so need to
		# select the derivatives functions to use. Could be simplified to single
		# functions that take a jit argument from solve_ivp in the future.
		if jit:
			derivatives = self.derivatives_jit
			derivatives_jacobian = self.derivatives_jacobian_jit
		else:
			derivatives = self.derivatives
			derivatives_jacobian = self.derivatives_jacobian

		# Note: odeint has issues solving with a long time step so need to use solve_ivp
		sol = integrate.solve_ivp(
			derivatives, [0, time_limit], y_init,
			method="LSODA", t_eval=[0, time_limit],
			jac=derivatives_jacobian)
		y = sol.y.T

		if np.any(y[-1, :] * (cellVolume * nAvogadro) <= -1):
			raise ValueError('Have negative values at equilibrium steady state -- probably due to numerical instability.')
		if np.linalg.norm(derivatives(0, y[-1, :]), np.inf) * (cellVolume * nAvogadro) > 1:
			raise RuntimeError('Did not reach steady state for equilibrium.')
		y[y < 0] = 0
		yMolecules = y * (cellVolume * nAvogadro)

		# Pick rounded solution that does not cause negative counts
		dYMolecules = yMolecules[-1, :] - yMolecules[0, :]
		for i in range(max_iter):
			rxnFluxes = stochasticRound(random_state, np.dot(self.mets_to_rxn_fluxes, dYMolecules))
			if np.all(moleculeCounts + self._stoichMatrix.dot(rxnFluxes) >= 0):
				break
		else:
			raise ValueError('Negative counts in equilibrium steady state.')

		rxnFluxesN = -1. * (rxnFluxes < 0) * rxnFluxes
		rxnFluxesP =  1. * (rxnFluxes > 0) * rxnFluxes
		moleculesNeeded = np.dot(self.Rp, rxnFluxesP) + np.dot(self.Pp, rxnFluxesN)

		return rxnFluxes, moleculesNeeded

	def get_monomers(self, cplxId):
		'''
		Returns subunits for a complex (or any ID passed). If the ID passed is
		already a monomer returns the monomer ID again with a stoichiometric
		coefficient of one.
		'''
		info = self._moleculeRecursiveSearch(cplxId, self._stoichMatrix, self.molecule_names)
		return {'subunitIds': np.array(list(info.keys())), 'subunitStoich': np.array(list(info.values()))}

	def get_metabolite(self, cplxId):
		D = self.get_monomers(cplxId)
		if len(D["subunitIds"]) > 2:
			raise Exception(
				"Calling this function only makes sense for reactions with 2 reactants")
		for subunit in D["subunitIds"]:
			if subunit in self.metabolite_set:
				return subunit

	def get_metabolite_coeff(self, cplxId):
		D = self.get_monomers(cplxId)
		if len(D["subunitIds"]) > 2:
			raise Exception(
				"Calling this function only makes sense for reactions with 2 reactants")
		for subunit, stoich in zip(D["subunitIds"], D["subunitStoich"]):
			if subunit in self.metabolite_set:
				return stoich

	def get_unbound(self, cplxId):
		D = self.get_monomers(cplxId)
		if len(D["subunitIds"]) > 2:
			raise Exception(
				"Calling this function only makes sense for reactions with 2 reactants")
		for subunit in D["subunitIds"]:
			if subunit not in self.metabolite_set:
				return subunit

	def get_fwd_rate(self, cplxId):
		rxnIdx = self.complex_name_to_rxn_idx[cplxId]
		return self.rates_fwd[rxnIdx]

	def get_rev_rate(self, cplxId):
		rxnIdx = self.complex_name_to_rxn_idx[cplxId]
		return self.rates_rev[rxnIdx]

	def set_fwd_rate(self, cplxId, rate):
		rxnIdx = self.complex_name_to_rxn_idx[cplxId]
		self.rates_fwd[rxnIdx] = rate

	def set_rev_rate(self, cplxId, rate):
		rxnIdx = self.complex_name_to_rxn_idx[cplxId]
		self.rates_rev[rxnIdx] = rate

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
			species = speciesList[i]

			if val != 0:
				x = self._moleculeRecursiveSearch(species, stoichMatrix, speciesList)
				for j in x:
					if j in total:
						total[j] += x[j]*(np.absolute(val))
					else:
						total[j] = x[j]*(np.absolute(val))
		return total
