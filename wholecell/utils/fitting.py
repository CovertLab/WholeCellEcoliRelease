from __future__ import absolute_import, division, print_function

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

def masses_and_counts_for_homeostatic_target(
		dry_mass_of_non_small_molecules,
		concentrations,
		weights,
		cell_density,
		avogadros_number
		):
	'''
	Computes the dry mass fractions and counts associated with small molecules to maintain
	concentrations consistent with targets.  (Also includes water.)

	The cell is composed of a number of 'mass fractions' i.e. DNA, RNA, protein, water, and
	the less specific "small molecules" which includes both inorganic and organic molecular species
	that form part of a cell.  While we take many of the former calculations as ground truth, we
	chose to adjust (recompute) the small molecule mass fraction according to per-molecule
	observations of small molecule concentrations (compiled from various sources).

	However, this creates a potential issue: we need the small molecule mass to compute the volume,
	and the volume in turn is used to compute the counts (and therefore masses) of the small
	molecules.  We denote the first small molecule mass as Ms, and the second as Ms'.

	The total mass of the cell, Mt, is the sum of the small and non-small molecule masses:

	Mt = Ms + Mns

	The volume of the cell V times the density of the cell rho is Mt, and therefore

	rho * V = Ms + Mns

	Ms = rho * V - Mns

	This gives us our first calculation of the small molecule mass.  For the second calculation, we
	first find the abundance of each small molecule species (count n_i) as

	n_i = V * c_i

	where c_i is the concentration of each species.  Then the mass associated with each species is

	m_i = V * w_i * c_i

	where w_i is the molecular weight of a given species,.  Finally, the total small molecule mass,
	estimated from small molecule counts, is

	Ms' = sum_i m_i = V * w^T c

	where w^T c is the dot-product between the two vectors.

	Equating Ms' and Ms, and solving for V:

	V = Mns / (rho - w^Tc)

	This allows us to compute the new volume, from which we can also compute and return all n_i and
	m_i.

	Parameters
	----------
	dry_mass_of_non_small_molecules: float unit'd scalar, dimensions of mass
		The total mass of the cell, minus the 'wet' mass (water) and the dry mass of other small
		molecules.
	concentrations: 1-D float unit'd array, with dimensions of concentration
		The target concentrations of the small molecules.
	weights: 1-D float unit'd array, with dimensions of mass per mol
		The molecular weights of the small molecules.
	cell_density: float unit'd scalar, dimensions of mass per volume
		The total density of the cell (wet and dry mass).
	avogadros_number: float unit'd scalar, dimensions of per mol
		The number of molecules per mole.

	Returns
	-------
	masses: 1-D float unit'd array, with dimensions of mass
		The mass associated with each molecular species.
	counts: 1-D float unit'd array, dimensionless (i.e. counts)
		The counts associated with each molecular species.
	'''

	# Compute the total mass concentration of the small molecules

	total_small_mol_mass_conc = np.dot(weights, concentrations)

	# Compute the new total cell volume that accomodates the small molecule concentrations

	cell_volume = dry_mass_of_non_small_molecules / (cell_density - total_small_mol_mass_conc)

	# Calculate and return the counts of molecules and their associated masses

	mols = cell_volume * concentrations

	counts = mols * avogadros_number
	masses = weights * mols

	return masses, counts

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

def cosine_similarity(samples):
	"""
	Finds the cosine similarity between samples.

	samples is a matrix of size (n_samples, sample_size)

	The output is a matrix of size (n_samples, n_samples).

	The cosine similarity is the normalized dot product between two
	vectors.  The name originates from the fact that the normalized dot
	product between two vectors is equal to the cosine of the angle
	formed by the two vectors.
	"""

	magnitudes = np.sqrt(np.sum(np.square(samples), 1))

	normed = samples / magnitudes[:, None]

	return normed.dot(normed.T)
