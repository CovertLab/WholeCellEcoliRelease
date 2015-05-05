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
	relativeExpression	- ndarray of floats	-	Relative expression of each species
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

def calcProteinCounts(kb, monomerMass):
	monomerExpression = calcProteinDistribution(kb)

	nMonomers = calcProteinTotalCounts(kb, monomerMass, monomerExpression)

	return nMonomers * monomerExpression


def calcProteinTotalCounts(kb, monomerMass, monomerExpression):
	return countsFromMassAndExpression(
		monomerMass.asNumber(units.g),
		kb.process.translation.monomerData["mw"].asNumber(units.g / units.mol),
		monomerExpression,
		kb.constants.nAvogadro.asNumber(1 / units.mol)
		)

def calcProteinDistribution(kb):
	return normalize(
		kb.process.transcription.rnaData['expression'][kb.relation.rnaIndexToMonomerMapping] /
		(np.log(2) / kb.doubling_time.asNumber(units.s) + kb.process.translation.monomerData["degRate"].asNumber(1 / units.s))
		)