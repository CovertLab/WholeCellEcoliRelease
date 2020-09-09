"""
Mass

Mass listener. Represents the total cellular mass.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

# TODO: generalize this logic for use with a generic simulation

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener
from wholecell.utils import units
import six

class Mass(wholecell.listeners.listener.Listener):
	""" Mass """

	_name = 'Mass'

	# Constructor
	def __init__(self, *args, **kwargs):
		# References to other internal states
		self.internal_states = None

		# NOTE: molecule weight is converted to femtograms/molecule from
		# grams/mol in BulkMolecules
		self.massUnits = 'fg'

		super(Mass, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Mass, self).initialize(sim, sim_data)

		self.internal_states = sim.internal_states

		self.processNames = list(sim.processes.keys())

		self.cellCycleLen = sim_data.condition_to_doubling_time[sim_data.condition].asNumber(units.s)

		self.rnaIndexes = np.array([
			sim_data.submass_name_to_index[name]
			for name in ["rRNA", "tRNA", "mRNA", "miscRNA", "nonspecific_RNA"]
			])

		self.rRnaIndex = sim_data.submass_name_to_index["rRNA"]
		self.smallMoleculeIndex = sim_data.submass_name_to_index["metabolite"]
		self.tRnaIndex = sim_data.submass_name_to_index["tRNA"]
		self.mRnaIndex = sim_data.submass_name_to_index["mRNA"]
		self.dnaIndex = sim_data.submass_name_to_index["DNA"]
		self.proteinIndex = sim_data.submass_name_to_index["protein"]
		self.waterIndex = sim_data.submass_name_to_index["water"]

		self.cellDensity = sim_data.constants.cell_density.asNumber(units.g / units.L)

		# Set initial values
		self.setInitial = False

		self.dryMass = 0.0
		# TODO: set initial masses based on some calculations of the expected
		# mother cell (divided by two) in the last time step

		# Register logged quantities
		self.registerLoggedQuantity(
			"Cell mass\n(fg)",
			"cellMass",
			".2f"
			)

		self.registerLoggedQuantity(
			"Dry mass\n(fg)",
			"dryMass",
			".2f"
			)

		self.registerLoggedQuantity(
			"Dry mass\nfold change",
			"dryMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"Expected\nfold change",
			"expectedMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"Growth\n(fg/s)",
			"growth",
			".4f"
			)

		self.registerLoggedQuantity(
			"Protein\nfraction",
			"proteinMassFraction",
			".3f"
			)

		self.registerLoggedQuantity(
			"Protein\nfold change",
			"proteinMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"RNA\nfraction",
			"rnaMassFraction",
			".3f"
			)

		self.registerLoggedQuantity(
			"RNA\nfold change",
			"rnaMassFoldChange",
			".3f"
			)

		self.registerLoggedQuantity(
			"Small mol\nfold change",
			"smallMoleculeFoldChange",
			".3f"
			)


	def update(self):
		oldDryMass = self.dryMass

		all_submasses = sum(
			state.mass() for state in six.viewvalues(self.internal_states))

		self.cellMass = all_submasses.sum()  # sum over all submasses

		self.waterMass = all_submasses[self.waterIndex]
		self.dryMass = self.cellMass - self.waterMass
		self.rnaMass = all_submasses[self.rnaIndexes].sum()
		self.rRnaMass = all_submasses[self.rRnaIndex]
		self.tRnaMass = all_submasses[self.tRnaIndex]
		self.mRnaMass = all_submasses[self.mRnaIndex]
		self.dnaMass = all_submasses[self.dnaIndex]
		self.proteinMass = all_submasses[self.proteinIndex]
		self.smallMoleculeMass = all_submasses[self.smallMoleculeIndex]

		# TODO (Eran) use this volume everywhere in the codebase that is currently calculating volume
		self.volume = self.cellMass / self.cellDensity

		self.processMassDifferences = sum(
			state.process_mass_diffs() for state in six.viewvalues(self.internal_states)
			).sum(axis=1)

		if self.simulationStep() > 0:
			self.growth = self.dryMass - oldDryMass

		else:
			self.growth = np.nan

		self.instantaniousGrowthRate = self.growth / self.timeStepSec() / self.dryMass

		self.proteinMassFraction = self.proteinMass / self.dryMass
		self.rnaMassFraction = self.rnaMass / self.dryMass

		if not self.setInitial:
			self.setInitial = True

			self.timeInitial = self.time()

			self.dryMassInitial = self.dryMass
			self.proteinMassInitial = self.proteinMass
			self.rnaMassInitial = self.rnaMass
			self.smallMoleculeMassInitial = self.smallMoleculeMass


		self.dryMassFoldChange = self.dryMass / self.dryMassInitial
		self.proteinMassFoldChange = self.proteinMass / self.proteinMassInitial
		self.rnaMassFoldChange = self.rnaMass / self.rnaMassInitial
		self.smallMoleculeFoldChange = self.smallMoleculeMass / self.smallMoleculeMassInitial

		self.expectedMassFoldChange = np.exp(np.log(2) * (self.time() - self.timeInitial) / self.cellCycleLen)


	def tableCreate(self, tableWriter):
		# Store units as metadata
		tableWriter.writeAttributes(
			cell_units = self.massUnits,
			cellDry_units = self.massUnits,
			metabolite_units = self.massUnits,
			growth_units = self.massUnits,
			rna_units = self.massUnits,
			protein_units = self.massUnits,
			water_units = self.massUnits,
			nucleoid_units = self.massUnits,
			processNames = self.processNames,
			smallMoleculeMass = self.smallMoleculeMass,
			subcolumns = {})


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			cellMass = self.cellMass,
			growth = self.growth,
			dryMass = self.dryMass,
			rnaMass = self.rnaMass,
			rRnaMass = self.rRnaMass,
			tRnaMass = self.tRnaMass,
			mRnaMass = self.mRnaMass,
			dnaMass = self.dnaMass,
			proteinMass = self.proteinMass,
			waterMass = self.waterMass,
			processMassDifferences = self.processMassDifferences.astype(np.float64),
			smallMoleculeMass = self.smallMoleculeMass,
			instantaniousGrowthRate = self.instantaniousGrowthRate,
			cellVolume = self.volume
			)
