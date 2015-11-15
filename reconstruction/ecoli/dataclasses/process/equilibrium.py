
from __future__ import division

import numpy as np
import os
import theano.tensor as T
import theano
import cPickle
import wholecell
from wholecell.utils import units
from . import metabolism
import scipy

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

		rxnIds = []

		stoichMatrixMass = []

		# Make sure reactions aren't duplicated in complexationReactions and equilibriumReactions
		equilibriumReactionIds = set([x["id"] for x in raw_data.equilibriumReactions])
		complexationReactionIds = set([x["id"] for x in raw_data.complexationReactions])

		if equilibriumReactionIds.intersection(complexationReactionIds) != set():
			raise Exception, "The following reaction ids are specified in equilibriumReactions and complexationReactions: %s" % (equilibriumReactionIds.intersection(complexationReactionIds))

		# Remove complexes that are currently not simulated
		FORBIDDEN_MOLECULES = {
			"modified-charged-selC-tRNA", # molecule does not exist
			}

		# Remove reactions that we know won't occur (e.g., don't do computations on metabolites that have zero counts)
		MOLECULES_THAT_WILL_EXIST_IN_SIMULATION = metabolism.METABOLITE_CONCENTRATIONS
		deleteReactions = []
		for reactionIndex, reaction in enumerate(raw_data.equilibriumReactions):
			for molecule in reaction["stoichiometry"]:
				if molecule["molecule"] in FORBIDDEN_MOLECULES or (molecule["type"] == "metabolite" and molecule["molecule"] not in MOLECULES_THAT_WILL_EXIST_IN_SIMULATION):
					deleteReactions.append(reactionIndex)
					break

		for reactionIndex in deleteReactions[::-1]:
			del raw_data.equilibriumReactions[reactionIndex]

		for reactionIndex, reaction in enumerate(raw_data.equilibriumReactions):
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

				if coefficient > 0:
					assert molecule["type"] == "proteincomplex"

				# Find molecular mass
				molecularMass = sim_data.getter.getMass([moleculeName]).asNumber(units.g / units.mol)[0]
				stoichMatrixMass.append(molecularMass)


		self._stoichMatrixI = np.array(stoichMatrixI)
		self._stoichMatrixJ = np.array(stoichMatrixJ)
		self._stoichMatrixV = np.array(stoichMatrixV)

		self.moleculeNames = molecules
		self.rxnIds = rxnIds
		self.ratesFwd = np.array(ratesFwd)
		self.ratesRev = np.array(ratesRev)

		# Mass balance matrix
		self._stoichMatrixMass = np.array(stoichMatrixMass)
		self.balanceMatrix = self.stoichMatrix()*self.massMatrix()

		# Find the mass balance of each equation in the balanceMatrix
		massBalanceArray = self.massBalance()

		# The stoichometric matrix should balance out to numerical zero.
		assert np.max([abs(x) for x in massBalanceArray]) < 1e-9

		self._makeMatrices()
		self._populateDerivativeAndJacobian()

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

	def _populateDerivativeAndJacobian(self):
		# TODO: Decide if this caching is worthwhile
		# TODO: Unhack this--this assumes a directory structure
		fixturesDir = os.path.join(
			os.path.dirname(os.path.dirname(wholecell.__file__)),
			"fixtures",
			"equilibrium"
			)

		needToCreate = False

		if not os.path.exists(fixturesDir):
			os.makedirs(fixturesDir)

		if os.path.exists(os.path.join(fixturesDir, "S.cPickle")):
			S = cPickle.load(open(os.path.join(fixturesDir, "S.cPickle"), "rb"))
			if not np.all(S == self.stoichMatrix()):
				needToCreate = True
		else:
			needToCreate = True

		if os.path.exists(os.path.join(fixturesDir, "ratesFwd.cPickle")):
			ratesFwd =  cPickle.load(open(os.path.join(fixturesDir, "ratesFwd.cPickle"), "rb"))
			if not np.all(ratesFwd == self.ratesFwd):
				needToCreate = True
		else:
			needToCreate = True

		if os.path.exists(os.path.join(fixturesDir, "ratesRev.cPickle")):
			ratesRev =  cPickle.load(open(os.path.join(fixturesDir, "ratesRev.cPickle"), "rb"))
			if not np.all(ratesRev == self.ratesRev):
				needToCreate = True
		else:
			needToCreate = True

		if needToCreate:
			self._makeMatrices()
			self._makeDerivative()
			cPickle.dump(self.stoichMatrix(), open(os.path.join(fixturesDir, "S.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(self.ratesFwd, open(os.path.join(fixturesDir, "ratesFwd.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(self.ratesRev, open(os.path.join(fixturesDir, "ratesRev.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(self.derivatives, open(os.path.join(fixturesDir, "derivatives.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(self.derivativesJacobian, open(os.path.join(fixturesDir, "jacobian.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
		else:
			self.derivatives = cPickle.load(open(os.path.join(fixturesDir, "derivatives.cPickle"), "rb"))
			self.derivativesJacobian = cPickle.load(open(os.path.join(fixturesDir, "jacobian.cPickle"), "rb"))

	def _makeMatrices(self):
		EPS = 1e-9

		S = self.stoichMatrix()
		S1 = np.zeros_like(S)
		S1[S < -1 * EPS] = -1
		S1[S > EPS] = 1

		Rp =  1. * (S1 < 0)
		Pp =  1. * (S1 > 0)

		self.Rp = Rp
		self.Pp = Pp

		metsToRxnFluxes = self.stoichMatrix().copy()

		metsToRxnFluxes[(np.abs(metsToRxnFluxes) > EPS).sum(axis = 1) > 1, : ] = 0
		for colIdx in xrange(metsToRxnFluxes.shape[1]):
			try:
				firstNonZeroIdx = np.where(np.abs(metsToRxnFluxes[:, colIdx]) > EPS)[0][0]
			except IndexError:
				raise Exception, "Column %d of S matrix not linearly independent!" % colIdx
			metsToRxnFluxes[:firstNonZeroIdx, colIdx] = 0
			metsToRxnFluxes[(firstNonZeroIdx + 1):, colIdx] = 0

		self.metsToRxnFluxes = metsToRxnFluxes.T

	def _makeDerivative(self):

		S = self.stoichMatrix()

		y = T.dvector()
		dy = [0 * y[0] for _ in xrange(S.shape[0])]
		for colIdx in xrange(S.shape[1]):
			negIdxs = np.where(S[:, colIdx] < 0)[0]
			posIdxs = np.where(S[:, colIdx] > 0)[0]

			reactantFlux = self.ratesFwd[colIdx]
			for negIdx in negIdxs:
				reactantFlux *= (y[negIdx] ** (-1 * S[negIdx, colIdx]))

			productFlux = self.ratesRev[colIdx]
			for posIdx in posIdxs:
				productFlux *=  (y[posIdx] ** ( 1 * S[posIdx, colIdx]))

			fluxForNegIdxs = (-1. * reactantFlux) + (1. * productFlux)
			fluxForPosIdxs = ( 1. * reactantFlux) - (1. * productFlux)

			for thisIdx in negIdxs:
				dy[thisIdx] += fluxForNegIdxs
			for thisIdx in posIdxs:
				dy[thisIdx] += fluxForPosIdxs

		t = T.dscalar()

		J = [T.grad(dy[i], y) for i in xrange(len(dy))]

		self.derivativesJacobian = theano.function([y, t], T.stack(*J), on_unused_input = "ignore")
		self.derivatives = theano.function([y, t], T.stack(*dy), on_unused_input = "ignore")


	# TODO: Should this method be here?
	# It could be useful in both the fitter and in the simulations
	# But it isn't just data
	def fluxesAndMoleculesToSS(self, moleculeCounts, cellVolume, nAvogadro):
		y_init = moleculeCounts / (cellVolume * nAvogadro)
		y = scipy.integrate.odeint(self.derivatives, y_init, t = [0, 1e4], Dfun = self.derivativesJacobian)

		if np.any(y[-1, :] * (cellVolume * nAvogadro) <= -1):
			raise Exception, "Have negative values -- probably due to numerical instability"
		if np.linalg.norm(self.derivatives(y[-1, :], 0), np.inf) * (cellVolume * nAvogadro) > 1:
			raise Exception, "Didn't reach steady state"
		y[y < 0] = 0
		yMolecules = y * (cellVolume * nAvogadro)

		dYMolecules = yMolecules[-1, :] - yMolecules[0, :]
		rxnFluxes = np.round(np.dot(self.metsToRxnFluxes, dYMolecules))
		rxnFluxesN = -1. * (rxnFluxes < 0) * rxnFluxes
		rxnFluxesP =  1. * (rxnFluxes > 0) * rxnFluxes
		moleculesNeeded = np.dot(self.Rp, rxnFluxesP) + np.dot(self.Pp, rxnFluxesN)
		return rxnFluxes, moleculesNeeded


	# TODO: These methods might not be necessary, consider deleting if that's the case
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
