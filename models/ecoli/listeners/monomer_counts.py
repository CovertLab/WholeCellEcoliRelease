#!/usr/bin/env python

"""
MonomerCounts Listener

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/15/2018
"""

from __future__ import division

import numpy as np
import wholecell.listeners.listener


class MonomerCounts(wholecell.listeners.listener.Listener):
	"""
	Listener for the full counts of protein monomers, including those that are
	part of a complex.
	"""
	_name = 'MonomerCounts'

	def __init__(self, *args, **kwargs):
		super(MonomerCounts, self).__init__(*args, **kwargs)

	def initialize(self, sim, sim_data):
		super(MonomerCounts, self).initialize(sim, sim_data)

		# Get IDs of all bulk molecules
		self.bulkMolecules = sim.internal_states["BulkMolecules"]
		bulk_molecule_ids = self.bulkMolecules.container.objectNames()

		# Get IDs of molecules involved in complexation and equilibrium
		complexation_molecule_ids = sim_data.process.complexation.moleculeNames
		complexation_complex_ids = sim_data.process.complexation.ids_complexes
		equilibrium_molecule_ids = sim_data.process.equilibrium.moleculeNames
		equilibrium_complex_ids = sim_data.process.equilibrium.ids_complexes
		self.monomer_ids = sim_data.process.translation.monomerData["id"].tolist()

		# Get IDs of ribosome subunits
		ribosome_50s_subunits = sim_data.process.complexation.getMonomers(
			sim_data.moleculeIds.s50_fullComplex)
		ribosome_30s_subunits = sim_data.process.complexation.getMonomers(
			sim_data.moleculeIds.s30_fullComplex)
		ribosome_subunit_ids = (ribosome_50s_subunits["subunitIds"].tolist() +
			ribosome_30s_subunits["subunitIds"].tolist())

		# Get IDs of RNA polymerase subunits
		rnap_subunits = sim_data.process.complexation.getMonomers(
			sim_data.moleculeIds.rnapFull)
		rnap_subunit_ids = rnap_subunits["subunitIds"].tolist()

		# Get IDs of replisome subunits
		replisome_trimer_subunits = sim_data.moleculeGroups.replisome_trimer_subunits
		replisome_monomer_subunits = sim_data.moleculeGroups.replisome_monomer_subunits
		replisome_subunit_ids = replisome_trimer_subunits + replisome_monomer_subunits

		# Get stoichiometric matrices for complexation, equilibrium, and the
		# assembly of unique molecules
		self.complexation_stoich = sim_data.process.complexation.stoichMatrixMonomers()
		self.equilibrium_stoich = sim_data.process.equilibrium.stoichMatrixMonomers()
		self.ribosome_stoich = np.hstack(
			(ribosome_50s_subunits["subunitStoich"],
			ribosome_30s_subunits["subunitStoich"]))
		self.rnap_stoich = rnap_subunits["subunitStoich"]
		self.replisome_stoich = np.hstack(
			(3*np.ones(len(replisome_trimer_subunits)),
			np.ones(len(replisome_monomer_subunits))))

		# Construct dictionary to quickly find bulk molecule indexes from IDs
		molecule_dict = {mol: i for i, mol in enumerate(bulk_molecule_ids)}

		def get_molecule_indexes(keys):
			return np.array([molecule_dict[x] for x in keys])

		# Get indexes of all relevant bulk molecules
		self.monomer_idx = get_molecule_indexes(self.monomer_ids)
		self.complexation_molecule_idx = get_molecule_indexes(complexation_molecule_ids)
		self.complexation_complex_idx = get_molecule_indexes(complexation_complex_ids)
		self.equilibrium_molecule_idx = get_molecule_indexes(equilibrium_molecule_ids)
		self.equilibrium_complex_idx = get_molecule_indexes(equilibrium_complex_ids)
		self.ribosome_subunit_idx = get_molecule_indexes(ribosome_subunit_ids)
		self.rnap_subunit_idx = get_molecule_indexes(rnap_subunit_ids)
		self.replisome_subunit_idx = get_molecule_indexes(replisome_subunit_ids)

		# Get indexes of all unique molecules that need to be accounted for
		self.uniqueMolecules = sim.internal_states["UniqueMolecules"]
		unique_molecule_ids = self.uniqueMolecules.container.objectNames()
		self.ribosome_idx = unique_molecule_ids.index('active_ribosome')
		self.rnap_idx = unique_molecule_ids.index('active_RNAP')
		self.replisome_idx = unique_molecule_ids.index("active_replisome")

	def allocate(self):
		super(MonomerCounts, self).allocate()

		self.monomerCounts = np.zeros(
			len(self.monomer_ids),
			np.int64
			)

	def update(self):
		# Get current counts of bulk and unique molecules
		bulkMoleculeCounts = self.bulkMolecules.container.counts()
		uniqueMoleculeCounts = self.uniqueMolecules.container.counts()
		n_active_ribosome = uniqueMoleculeCounts[self.ribosome_idx]
		n_active_rnap = uniqueMoleculeCounts[self.rnap_idx]
		n_active_replisome = uniqueMoleculeCounts[self.replisome_idx]

		# Account for monomers in bulk molecule complexes
		complex_monomer_counts = np.dot(self.complexation_stoich,
			np.negative(bulkMoleculeCounts[self.complexation_complex_idx]))
		equilibrium_monomer_counts = np.dot(self.equilibrium_stoich,
			np.negative(bulkMoleculeCounts[self.equilibrium_complex_idx]))

		bulkMoleculeCounts[self.complexation_molecule_idx] += complex_monomer_counts.astype(np.int)
		bulkMoleculeCounts[self.equilibrium_molecule_idx] += equilibrium_monomer_counts.astype(np.int)

		# Account for monomers in unique molecule complexes
		n_ribosome_subunit = n_active_ribosome * self.ribosome_stoich
		n_rnap_subunit = n_active_rnap * self.rnap_stoich
		n_replisome_subunit = n_active_replisome * self.replisome_stoich
		bulkMoleculeCounts[self.ribosome_subunit_idx] += n_ribosome_subunit.astype(np.int)
		bulkMoleculeCounts[self.rnap_subunit_idx] += n_rnap_subunit.astype(np.int)
		bulkMoleculeCounts[self.replisome_subunit_idx] += n_replisome_subunit.astype(np.int)

		# Update monomerCounts
		self.monomerCounts = bulkMoleculeCounts[self.monomer_idx]

	def tableCreate(self, tableWriter):
		subcolumns = {
			'monomerCounts': 'monomerIds'}

		tableWriter.writeAttributes(
			monomerIds = self.monomer_ids,
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			monomerCounts = self.monomerCounts,
			)
