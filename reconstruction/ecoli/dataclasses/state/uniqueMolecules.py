"""
SimulationData for unique molecules state

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/19/2015
"""

from __future__ import division
import collections

import numpy as np

from wholecell.utils import units
from wholecell.utils.unit_struct_array import UnitStructArray

class UniqueMolecules(object):
	""" UniqueMolecules """

	def __init__(self, raw_data, sim_data):
		self._buildUniqueMolecules(raw_data, sim_data)


	def _buildUniqueMolecules(self, raw_data, sim_data):
		self.uniqueMoleculeDefinitions = collections.OrderedDict([
			("activeRnaPoly", {
				'rnaIndex' : 'i8',
				'transcriptLength' : 'i8'
				}),
			("activeRibosome", {
				'proteinIndex' : 'i8',
				'peptideLength': 'i8'
				}),

			("dnaPolymerase", {
				'chromosomeLocation' : 'i8',
				'directionIsPositive' : 'bool',
				'isLeading' : 'bool'
				}),
			])

		rnaPolyComplexMass = self.bulkMolecules["mass"][self.bulkMolecules["moleculeId"] == "APORNAP-CPLX[c]"].asNumber()

		# TODO: This is a bad hack that works because in the fitter
		# I have forced expression to be these subunits only
		ribosomeSubunits = [self.s30_fullComplex, self.s50_fullComplex]

		ribosomeMass = sum(
			entry["mass"] for entry in self.bulkMolecules.struct_array
			if entry["moleculeId"] in ribosomeSubunits
			)

		dnaPolyMass = np.zeros_like(rnaPolyComplexMass) # NOTE: dnaPolymerases currently have no mass

		uniqueMoleculeMasses = np.zeros(
			shape = len(self.uniqueMoleculeDefinitions),
			dtype = [
				('moleculeId', 'a50'),
				("mass", "{}f8".format(len(MOLECULAR_WEIGHT_ORDER))),
				]
			)

		uniqueMoleculeMasses["moleculeId"] = self.uniqueMoleculeDefinitions.keys()

		uniqueMoleculeMasses["mass"] = np.vstack([
			rnaPolyComplexMass,
			ribosomeMass,
			dnaPolyMass
			])

		self.uniqueMoleculeMasses = UnitStructArray(
			uniqueMoleculeMasses,
			{"moleculeId":None, "mass":units.g / units.mol}
			)

		# TODO: add the ability to "register" a bulk molecule as a unique
		# molecule to handle most of the above logic