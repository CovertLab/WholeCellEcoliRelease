"""
Two component systems.

Note: Ligand binding to histidine kinases is modeled by equilibrium.

TODOs:
_populateDerivativeAndJacobian()
	Decide if this caching is worthwhile
	Assumes a directory structure

moleculesToNextTimeStep()
	Consider relocating (since it's useful for both the parca and simulation)

"""

from __future__ import absolute_import, division, print_function

import numpy as np
import os
import cPickle
import scipy
import re
import sympy as sp

import wholecell
from wholecell.utils import filepath
from wholecell.utils import units
from wholecell.utils.write_ode_file import writeOdeFile


class TwoComponentSystem(object):
	def __init__(self, raw_data, sim_data):
		# Store two component system raw data for use in analysis
		sim_data.moleculeGroups.twoComponentSystems = raw_data.twoComponentSystems

		# Build the abstractions needed for two component systems
		molecules = []
		moleculeTypes = []

		ratesFwd = []
		ratesRev = []
		rxnIds = []

		stoichMatrixI = []
		stoichMatrixJ = []
		stoichMatrixV = []

		stoichMatrixMass = []

		independentMolecules = []
		independent_molecule_indexes = []
		independentToDependentMolecules = {}

		activeToInactiveTF = {} #convention: active TF is the DNA-binding form

		# Build template reactions
		signalingTemplate = {
			1: ["POS-LIGAND-BOUND-HK-PHOSPHORYLATION_RXN",
				"POS-LIGAND-BOUND-HK-PHOSPHOTRANSFER_RXN",
				"POS-RR-DEPHOSPHORYLATION_RXN",
				"POS-HK-PHOSPHORYLATION_RXN",
				"POS-HK-PHOSPHOTRANSFER_RXN",
				],
			-1: ["NEG-HK-PHOSPHORYLATION_RXN",
				"NEG-HK-PHOSPHOTRANSFER_RXN",
				"NEG-RR-DEPHOSPHORYLATION_RXN",
				],
			}

		reactionTemplate = {}
		for reactionIndex, reaction in enumerate(raw_data.twoComponentSystemTemplates):
			reactionTemplate[str(reaction["id"])] = reaction

		# Build stoichiometry matrix
		for systemIndex, system in enumerate(raw_data.twoComponentSystems):
			for reaction in signalingTemplate[system["orientation"]]:
				reactionName = self.getReactionName(reaction, system["molecules"])

				if reactionName not in rxnIds:
					rxnIds.append(reactionName)
					ratesFwd.append(reactionTemplate[reaction]["forward rate"])
					ratesRev.append(reactionTemplate[reaction]["reverse rate"])
					reactionIndex = len(rxnIds) - 1

				else:
					reactionIndex = rxnIds.index(reactionName)

				for molecule in reactionTemplate[reaction]["stoichiometry"]:

					# Build name for system molecules
					if molecule["molecule"] in system["molecules"]:
						moleculeName = "{}[{}]".format(
							system["molecules"][molecule["molecule"]],
							molecule["location"]
							)

					# Build name for common molecules (ATP, ADP, PI, WATER, PROTON)
					else:
						moleculeName = "{}[{}]".format(
							molecule["molecule"],
							molecule["location"]
							)

					if moleculeName not in molecules:
						molecules.append(moleculeName)
						moleculeIndex = len(molecules) - 1
						moleculeTypes.append(str(molecule["molecule"]))

					else:
						moleculeIndex = molecules.index(moleculeName)

					coefficient = molecule["coeff"]

					# Store indices for the row and column, and molecule coefficient for building the stoichiometry matrix
					stoichMatrixI.append(moleculeIndex)
					stoichMatrixJ.append(reactionIndex)
					stoichMatrixV.append(coefficient)

					# Build matrix with linearly independent rows based on network orientation
					if str(molecule["molecule"]) in ["HK", "RR", "ATP"] and moleculeName not in independentMolecules:
						independentMolecules.append(moleculeName)
						independent_molecule_indexes.append(moleculeIndex)

						# Map linearly independent molecules (rows) to their dependents (phosphorylated forms of histidine kinases and response regulators)
						if str(molecule["molecule"]) != "ATP":
							independentToDependentMolecules[moleculeName] = "{}[{}]".format(
								system["molecules"]["PHOSPHO-" + str(molecule["molecule"])],
								molecule["location"]
								)

						# Map active transcription factors (phosphorylated response regulators) to their inactive forms (unphosphorylated response regulators)
						if str(molecule["molecule"]) == "RR":
							activeTF = "{}[{}]".format(
								system["molecules"]["PHOSPHO-RR"],
								molecule["location"]
								)
							activeToInactiveTF[activeTF] = moleculeName

					# Account for ligand-bound histidine kinases for positively oriented networks
					if system["orientation"] == 1:
						if str(molecule["molecule"]) == "HK-LIGAND" and moleculeName not in independentMolecules:
							independentMolecules.append(moleculeName)
							independent_molecule_indexes.append(moleculeIndex)

							# Map the linearly independent ligand-bound histidine kinases to their dependents (phosphorylated forms of ligand-bound histidine kinases)
							independentToDependentMolecules[moleculeName] = "{}[{}]".format(
								system["molecules"]["PHOSPHO-" + str(molecule["molecule"])],
								molecule["location"]
								)

					# Find molecular mass
					molecularMass = sim_data.getter.getMass([moleculeName]).asNumber(units.g / units.mol)[0]
					stoichMatrixMass.append(molecularMass)

		self._stoichMatrixI = np.array(stoichMatrixI)
		self._stoichMatrixJ = np.array(stoichMatrixJ)
		self._stoichMatrixV = np.array(stoichMatrixV)

		self.moleculeNames = np.array(molecules)
		self.moleculeTypes = np.array(moleculeTypes)
		self.rxnIds = rxnIds
		self.ratesFwd = np.array(ratesFwd)
		self.ratesRev = np.array(ratesRev)

		self.independentMolecules = np.array(independentMolecules)
		self.independent_molecule_indexes = np.array(independent_molecule_indexes)
		self.independentToDependentMolecules = independentToDependentMolecules

		self.independentMoleculesAtpIndex = np.where(self.independentMolecules == "ATP[c]")[0][0]

		self.complexToMonomer = self._buildComplexToMonomer(raw_data.modifiedFormsStoichiometry, self.moleculeNames)

		# Mass balance matrix
		self._stoichMatrixMass = np.array(stoichMatrixMass)
		self.balanceMatrix = self.stoichMatrix()*self.massMatrix()

		# Find the mass balance of each equation in the balanceMatrix
		massBalanceArray = self.massBalance()

		# The stoichometric matrix should balance out to numerical zero.
		assert np.max([abs(x) for x in massBalanceArray]) < 1e-9

		# Build matrices
		self._populateDerivativeAndJacobian()
		self.dependencyMatrix = self.makeDependencyMatrix()

		# Map active TF to inactive TF
		self.activeToInactiveTF = activeToInactiveTF


	def _buildComplexToMonomer(self, modifiedFormsMonomers, tcsMolecules):
		'''
		Maps each complex to a dictionary that maps each subunit of the complex to its stoichiometry
		'''
		D = {}
		for row in modifiedFormsMonomers:
			if str(row["complexID"]) in tcsMolecules:
				D[str(row["complexID"])] = {}
				for subunit in row["subunits"]:
					D[str(row["complexID"])][str(subunit["monomer"])] = float(subunit["stoichiometry"])

		return D


	def stoichMatrix(self):
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


	def massMatrix(self):
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


	def massBalance(self):
		'''
		Sum along the columns of the massBalance matrix to check for reaction mass balance
		'''
		return np.sum(self.balanceMatrix, axis=0)


	def stoichMatrixMonomers(self):
		'''
		Builds stoichiometry matrix for monomers (complex subunits)
		Rows: molecules (complexes and monomers) 
		Columns: complexes 
		Values: monomer stoichiometry
		'''
		ids_complexes = self.complexToMonomer.keys()
		stoichMatrixMonomersI = []
		stoichMatrixMonomersJ = []
		stoichMatrixMonomersV = []
		for colIdx, id_complex in enumerate(ids_complexes):
			D = self.getMonomers(id_complex)

			rowIdx = self.moleculeNames.tolist().index(id_complex)
			stoichMatrixMonomersI.append(rowIdx)
			stoichMatrixMonomersJ.append(colIdx)
			stoichMatrixMonomersV.append(1.)

			for subunitId, subunitStoich in zip(D["subunitIds"], D["subunitStoich"]):
				if subunitId in self.moleculeNames.tolist():		
					rowIdx = self.moleculeNames.tolist().index(subunitId)
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
		'''
		Creates callable functions for computing the derivative and the
		Jacobian.
		'''
		fixturesDir = filepath.makedirs(
			os.path.dirname(os.path.dirname(wholecell.__file__)),
			"fixtures",
			"twoComponentSystem"
			)
		odeFile = os.path.join(
			os.path.dirname(os.path.dirname(wholecell.__file__)),
			"reconstruction", "ecoli", "dataclasses", "process", "two_component_system_odes.py"
			)
		odeParcaFile = os.path.join(
			os.path.dirname(os.path.dirname(wholecell.__file__)),
			"reconstruction", "ecoli", "dataclasses", "process", "two_component_system_odes_parca.py"
			)

		needToCreate = False

		if not os.path.exists(odeFile):
			needToCreate = True

		if not os.path.exists(odeParcaFile):
			needToCreate = True

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
			self._makeDerivative()
			self._makeDerivativeParca()

			writeOdeFile(odeFile, self.derivativesSymbolic, self.derivativesJacobianSymbolic)
			writeOdeFile(odeParcaFile, self.derivativesParcaSymbolic, self.derivativesParcaJacobianSymbolic)

			# Modules are imported here to ensure the files exist before import
			import reconstruction.ecoli.dataclasses.process.two_component_system_odes
			import reconstruction.ecoli.dataclasses.process.two_component_system_odes_parca

			self.derivatives = reconstruction.ecoli.dataclasses.process.two_component_system_odes.derivatives
			self.derivatives_jacobian = reconstruction.ecoli.dataclasses.process.two_component_system_odes.derivativesJacobian
			self.derivatives_parca = reconstruction.ecoli.dataclasses.process.two_component_system_odes_parca.derivatives
			self.derivatives_parca_jacobian = reconstruction.ecoli.dataclasses.process.two_component_system_odes_parca.derivativesJacobian

			cPickle.dump(self.stoichMatrix(), open(os.path.join(fixturesDir, "S.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(self.ratesFwd, open(os.path.join(fixturesDir, "ratesFwd.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
			cPickle.dump(self.ratesRev, open(os.path.join(fixturesDir, "ratesRev.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
		else:
			# Modules are imported here to ensure the files exist before import
			import reconstruction.ecoli.dataclasses.process.two_component_system_odes
			import reconstruction.ecoli.dataclasses.process.two_component_system_odes_parca

			self.derivatives = reconstruction.ecoli.dataclasses.process.two_component_system_odes.derivatives
			self.derivatives_jacobian = reconstruction.ecoli.dataclasses.process.two_component_system_odes.derivativesJacobian
			self.derivatives_parca = reconstruction.ecoli.dataclasses.process.two_component_system_odes_parca.derivatives
			self.derivatives_parca_jacobian = reconstruction.ecoli.dataclasses.process.two_component_system_odes_parca.derivativesJacobian


	def _makeDerivative(self):
		'''
		Creates symbolic representation of the ordinary differential equations
		and the Jacobian. Used during simulations.
		'''
		S = self.stoichMatrix()

		yStrings = ["y[%d]" % x for x in xrange(S.shape[0])]
		y = sp.symbols(yStrings)
		dy = [sp.symbol.S.Zero] * S.shape[0]

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

		dy = sp.Matrix(dy)
		J = dy.jacobian(y)

		self.derivativesJacobianSymbolic = J
		self.derivativesSymbolic = dy


	def _makeDerivativeParca(self):
		'''
		Creates symbolic representation of the ordinary differential equations
		and the Jacobian assuming ATP, ADP, Pi, water and protons are at
		steady state. Used in the parca.
		'''
		S = self.stoichMatrix()

		yStrings = ["y[%d]" % x for x in xrange(S.shape[0])]
		y = sp.symbols(yStrings)
		dy = [sp.symbol.S.Zero] * S.shape[0]

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

		# Metabolism will keep these molecules at steady state
		constantMolecules = ["ATP[c]", "ADP[c]", "PI[c]", "WATER[c]", "PROTON[c]"]
		for molecule in constantMolecules:
			moleculeIdx = np.where(self.moleculeNames == molecule)[0][0]
			dy[moleculeIdx] = sp.symbol.S.Zero

		dy = sp.Matrix(dy)
		J = dy.jacobian(y)

		self.derivativesParcaJacobianSymbolic = J
		self.derivativesParcaSymbolic = dy


	def moleculesToNextTimeStep(self, moleculeCounts, cellVolume,
			nAvogadro, timeStepSec, solver="LSODA"):
		"""
		Calculates the changes in the counts of molecules in the next timestep
		by solving an initial value ODE problem.

		Args:
			moleculeCounts (1d ndarray, ints): current counts of molecules
			involved in the ODE
			cellVolume (float): current volume of cell
			nAvogadro (float): Avogadro's number
			timeStepSec (float): current length of timestep in seconds
			solver (str): name of the ODE solver to use

		Returns:
			moleculesNeeded (1d ndarray, ints): counts of molecules that need
			to be consumed
			allMoleculesChanges (1d ndarray, ints): expected changes in
			molecule counts after timestep
		"""
		y_init = moleculeCounts / (cellVolume * nAvogadro)

		if solver == "BDF":
			# Note: solve_ivp requires the order of arguments (t and y) for the
			# derivative and jacobian functions to be flipped relative to the
			# requirements of odeint. Wrapper functions were used to do this
			# without changing the original functions.
			sol = scipy.integrate.solve_ivp(
				self.derivatives_flipped, [0, timeStepSec], y_init,
				method="BDF", t_eval=[0, timeStepSec],
				jac=self.derivatives_jacobian_flipped
				)
			y = sol.y.T
		else:
			y = scipy.integrate.odeint(
				self.derivatives, y_init,
				t=[0, timeStepSec], Dfun=self.derivatives_jacobian
				)

		if np.any(y[-1, :] * (cellVolume * nAvogadro) <= -1):
			raise Exception(
				"Solution to ODE for two-component systems has negative values."
				)

		y[y < 0] = 0
		yMolecules = y * (cellVolume * nAvogadro)
		dYMolecules = yMolecules[-1, :] - yMolecules[0, :]

		independentMoleculesCounts = np.array(
			[np.round(dYMolecules[x]) for x in self.independent_molecule_indexes]
			)

		# To ensure that we have non-negative counts of phosphate, we must
		# have the following (which can be seen from the dependency matrix)
		independentMoleculesCounts[self.independentMoleculesAtpIndex] = (
			independentMoleculesCounts[:self.independentMoleculesAtpIndex].sum()
			+ independentMoleculesCounts[(self.independentMoleculesAtpIndex + 1):].sum()
			)

		# Calculate changes in molecule counts for all molecules
		allMoleculesChanges = np.dot(
			self.dependencyMatrix, independentMoleculesCounts)

		moleculesNeeded = np.negative(allMoleculesChanges).clip(min=0)

		return moleculesNeeded, allMoleculesChanges


	def moleculesToSS(self, moleculeCounts, cellVolume, nAvogadro, timeStepSec):
		"""
		Calculates the changes in the counts of molecules as the system
		reaches steady state

		Args:
			moleculeCounts: current counts of molecules involved in the ODE
			cellVolume: current volume of cell
			nAvogadro: Avogadro's number
			timeStepSec: current length of timestep (set to large number)

		Returns:
			moleculesNeeded: counts of molecules that need to be consumed
			allMoleculesChanges: expected changes in molecule counts after
			timestep
		"""
		# TODO (Gwanggyu): This function should probably get merged with the
		# function above.

		y_init = moleculeCounts / (cellVolume * nAvogadro)

		y = scipy.integrate.odeint(
			self.derivatives_parca, y_init,
			t=[0, timeStepSec], Dfun=self.derivatives_parca_jacobian
			)

		if np.any(y[-1, :] * (cellVolume * nAvogadro) <= -1):
			raise Exception(
				"Solution to ODE for two-component systems has negative values."
				)

		y[y < 0] = 0
		yMolecules = y * (cellVolume * nAvogadro)
		dYMolecules = yMolecules[-1, :] - yMolecules[0, :]

		independentMoleculesCounts = np.array(
			[np.round(dYMolecules[x]) for x in self.independent_molecule_indexes]
			)

		# To ensure that we have non-negative counts of phosphate, we must
		# have the following (which can be seen from the dependency matrix)
		independentMoleculesCounts[self.independentMoleculesAtpIndex] = (
			independentMoleculesCounts[:self.independentMoleculesAtpIndex].sum()
			+ independentMoleculesCounts[(self.independentMoleculesAtpIndex + 1):].sum()
			)

		# Calculate changes in molecule counts for all molecules
		allMoleculesChanges = np.dot(
			self.dependencyMatrix, independentMoleculesCounts)

		moleculesNeeded = np.negative(allMoleculesChanges).clip(min=0)

		return moleculesNeeded, allMoleculesChanges


	def getMonomers(self, cplxId):
		'''
		Returns subunits for a complex (or any ID passed).
		If the ID passed is already a monomer returns the
		monomer ID again with a stoichiometric coefficient
		of zero.
		'''

		info = self.complexToMonomer
		if cplxId in info:
			out = {'subunitIds' : info[cplxId].keys(), 'subunitStoich' : info[cplxId].values()}
		else:
			out = {'subunitIds' : cplxId, 'subunitStoich' : 1}
		return out


	def getReactionName(self, templateName, systemMolecules):
		'''
		Returns reaction name for a particular system.
		'''
		startIndex = 0
		reactionName = templateName
		for endIndex in [x.start() for x in re.finditer("-", templateName)]:
			if templateName[startIndex:endIndex] in systemMolecules:
				reactionName = reactionName.replace(templateName[startIndex:endIndex], str(systemMolecules[templateName[startIndex:endIndex]]))
			startIndex = endIndex + 1

		return reactionName


	def makeDependencyMatrix(self):
		'''
		Builds matrix mapping linearly independent molecules (ATP, histidine kinases, 
		response regulators, and ligand-bound histidine kinases for positively oriented 
		networks) to their dependents.
		'''
		moleculeTypes = self.moleculeTypes
		dependencyMatrixI = []
		dependencyMatrixJ = []
		dependencyMatrixV = []

		for independentMoleculeIndex, independentMoleculeId in enumerate(self.independent_molecule_indexes):
			dependencyMatrixI.append(independentMoleculeId)
			dependencyMatrixJ.append(independentMoleculeIndex)
			dependencyMatrixV.append(1)

			if self.moleculeNames[independentMoleculeId] == "ATP[c]":
				dependencyMatrixATPJ = independentMoleculeIndex
			else:
				moleculeType = moleculeTypes[independentMoleculeId]
				dependentMoleculeId = int(np.where(self.moleculeNames == self.independentToDependentMolecules[self.moleculeNames[independentMoleculeId]])[0])
				dependencyMatrixI.append(dependentMoleculeId)
				dependencyMatrixJ.append(independentMoleculeIndex)
				dependencyMatrixV.append(-1)

		# ATP dependents: ADP, PI, WATER, PROTON)
		for ATPdependent in ["ADP[c]", "PI[c]", "WATER[c]", "PROTON[c]"]:
			dependencyMatrixI.append(int(np.where(self.moleculeNames == ATPdependent)[0]))
			dependencyMatrixJ.append(dependencyMatrixATPJ)
			if ATPdependent == "WATER[c]":
				dependencyMatrixV.append(1)
			else:
				dependencyMatrixV.append(-1)

		for col in np.arange(self.independent_molecule_indexes.size):
			if col == dependencyMatrixATPJ:
				continue
			else:
				dependencyMatrixI.append(int(np.where(self.moleculeNames == "PI[c]")[0]))
				dependencyMatrixJ.append(col)
				dependencyMatrixV.append(1)

				dependencyMatrixI.append(int(np.where(self.moleculeNames == "WATER[c]")[0]))
				dependencyMatrixJ.append(col)
				dependencyMatrixV.append(-1)

		shape = (np.max(dependencyMatrixI) +1, np.max(dependencyMatrixJ) +1)

		out = np.zeros(shape, np.float64)

		out[dependencyMatrixI, dependencyMatrixJ] = dependencyMatrixV
		return out


	def derivatives_flipped(self, t, y):
		"""
		Wrapper function to flip the order of arguments of the derivative
		function given as an argument to solve_ivp.
		"""
		return self.derivatives(y, t)


	def derivatives_jacobian_flipped(self, t, y):
		"""
		Wrapper function to flip the order of arguments of the Jacobian
		function given as an argument to solve_ivp.
		"""
		return self.derivatives_jacobian(y, t)
