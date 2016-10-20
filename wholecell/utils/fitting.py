import numpy as np
import unum # Imported here to be used in getCountsFromMassAndExpression assertions
from wholecell.utils import units

def normalize(array):
	return np.array(array).astype("float") / np.linalg.norm(array, 1)

def countsFromMassAndExpression(mass, mws, relativeExpression, nAvogadro):
	"""
	countsFromMassAndExpression

	mass 				- float -				Total mass you want counts to sum to
	mws					- ndarray of floats -	Molecular weights of each species
	relativeExpression	- ndarray of floats	-	Relative expres`sion of each species
	nAvogadro 			- float -				Avogadro's number

	Example:
		mass = 10.
		mws = [10., 5.]
		relativeExpression = [0.33, 0.66]
		countsFromMassAndExpression(mass, mws, relativeExpression, nAvogadro) = 1.93e23
	"""
	assert np.allclose(np.sum(relativeExpression), 1)
	assert type(mass) != unum.Unum
	assert type(mws) != unum.Unum
	assert type(relativeExpression) != unum.Unum
	assert type(nAvogadro) != unum.Unum
	return mass / np.dot(mws / nAvogadro, relativeExpression)

def massesAndCountsToAddForHomeostaticTargets(massInitial, poolIds, poolConcentrations, mws, cellDensity, nAvogadro):
	diag = (cellDensity / (mws * poolConcentrations) - 1).asNumber()
	A = -1 * np.ones((diag.size, diag.size))
	A[np.diag_indices(diag.size)] = diag
	b = massInitial.asNumber(units.g) * np.ones(diag.size)

	massesToAdd = units.g * np.linalg.solve(A, b)
	countsToAdd = massesToAdd / mws * nAvogadro

	V = (massInitial + units.sum(massesToAdd)) / cellDensity

	assert np.allclose(
		(countsToAdd / nAvogadro / V).asNumber(units.mol / units.L),
		(poolConcentrations).asNumber(units.mol / units.L)
		)

	return massesToAdd, countsToAdd

def calcProteinCounts(sim_data, monomerMass):
	monomerExpression = calcProteinDistribution(sim_data)

	nMonomers = calcProteinTotalCounts(sim_data, monomerMass, monomerExpression)

	return nMonomers * monomerExpression


def calcProteinTotalCounts(sim_data, monomerMass, monomerExpression):
	return countsFromMassAndExpression(
		monomerMass.asNumber(units.g),
		sim_data.process.translation.monomerData["mw"].asNumber(units.g / units.mol),
		monomerExpression,
		sim_data.constants.nAvogadro.asNumber(1 / units.mol)
		)

def calcProteinDistribution(sim_data):
	return normalize(
		sim_data.process.transcription.rnaData["expression"][sim_data.relation.rnaIndexToMonomerMapping] /
		(np.log(2) / sim_data.doubling_time.asNumber(units.s) + sim_data.process.translation.monomerData["degRate"].asNumber(1 / units.s))
		)