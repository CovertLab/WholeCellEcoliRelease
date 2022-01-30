"""
Dynamics of limitations related to growth.  Excess ribosomal RNA or proteins
compared to fully formed ribosomes can show if more RNA or more protein needs
to be produced.  Comparing the ribosome fraction to amino acid synthesis fraction
can also show if supply or demand is limiting.

TODO:
- move function to another file (more general analysis helpers)?
- single analysis or something else?
- rename file if enzyme fraction doesn't make sense here
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import singleAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_bulk_molecules, read_stacked_columns
from wholecell.io.tablereader import TableReader
from wholecell.utils import units


def calculate_ribosome_excesses(sim_data, paths):
	# TODO: separate out sim_data processing from data loading for plots that might iterate over multiple variants/other groups
	complexation = sim_data.process.complexation
	metabolism = sim_data.process.metabolism
	rna_data = sim_data.process.transcription.rna_data

	# Get complexation data
	complex_stoich = -complexation.stoich_matrix_monomers().astype(int)
	monomer_ids = {id_: i for i, id_ in enumerate(complexation.molecule_names)}
	complex_ids = {id_: i for i, id_ in enumerate(complexation.ids_complexes)}

	# Get molecules of interest
	rnas = [rna for rna in rna_data['id'][rna_data['is_rRNA']] if rna in monomer_ids]
	proteins = sim_data.molecule_groups.ribosomal_proteins
	complexes = ([sim_data.molecule_ids.s50_full_complex]
		+ [sim_data.molecule_ids.s30_full_complex]
		+ sim_data.molecule_groups.s50_protein_complexes)
	ribosome_subunit_indices = 2  # the indices in complexes that directly form the ribosome (ie 30S and 50S)

	# Enzymes involved in mechanistic amino acid synthesis
	synthesis_enzymes = metabolism.aa_enzymes[metabolism.enzyme_to_amino_acid_fwd.sum(1).astype(bool)]
	synthesis_monomers = sorted({
		subunit
		for enzyme in synthesis_enzymes
		for subunit in complexation.get_monomers(enzyme)['subunitIds']
		})

	# Select subset of the complexation stoich matrix
	rna_idx = np.array([monomer_ids[id_] for id_ in rnas])
	protein_idx = np.array([monomer_ids[id_] for id_ in proteins])
	complex_idx = np.array([complex_ids[id_] for id_ in complexes])
	rna_stoich = complex_stoich[rna_idx, :][:, complex_idx].T
	protein_stoich = complex_stoich[protein_idx, :][:, complex_idx].T

	# Molecular weights for molecule groups
	mw_rnas = sim_data.getter.get_masses(rnas).asNumber(units.fg / units.count)
	mw_proteins = sim_data.getter.get_masses(proteins).asNumber(units.fg / units.count)
	mw_ribosome = (sim_data.getter.get_mass(sim_data.molecule_ids.s50_full_complex) + sim_data.getter.get_mass(sim_data.molecule_ids.s30_full_complex)).asNumber(units.fg / units.count)
	mw_enzymes = sim_data.getter.get_masses(synthesis_monomers).asNumber(units.fg / units.count)

	ribosome_mass_rna = (rna_stoich @ mw_rnas)[1:].sum()
	ribosome_mass_protein = (protein_stoich @ mw_proteins)[1:].sum()
	ribosome_fraction_rna = ribosome_mass_rna / (ribosome_mass_rna + ribosome_mass_protein)
	ribosome_fraction_protein = ribosome_mass_protein / (ribosome_mass_rna + ribosome_mass_protein)

	unique_molecule_reader = TableReader(os.path.join(paths[0], 'simOut', 'UniqueMoleculeCounts'))
	unique_molecule_ids = unique_molecule_reader.readAttribute('uniqueMoleculeIds')
	ribosome_idx = unique_molecule_ids.index('active_ribosome')

	monomers_reader = TableReader(os.path.join(os.path.join(paths[0], 'simOut', 'MonomerCounts')))
	protein_ids = {monomer: i for i, monomer in enumerate(monomers_reader.readAttribute('monomerIds'))}
	synthesis_idx = np.array([protein_ids[m] for m in synthesis_monomers])

	total_protein_mass = read_stacked_columns(paths, 'Mass', 'proteinMass',
		remove_first=True).squeeze()
	enzyme_mass = read_stacked_columns(paths, 'MonomerCounts', 'monomerCounts',
		remove_first=True, fun=lambda x: (x[:, synthesis_idx] @ mw_enzymes).reshape(-1, 1)).squeeze()
	active_ribosome_counts = read_stacked_columns(paths, 'UniqueMoleculeCounts', 'uniqueMoleculeCounts',
		remove_first=True, fun=lambda x: x[:, ribosome_idx].reshape(-1, 1)).squeeze()
	rna_counts, protein_counts, complex_counts = read_stacked_bulk_molecules(
		paths, (rnas, proteins, complexes), remove_first=True)

	# Add complex counts to RNA/protein
	inactive_ribosomes = complex_counts[:, :ribosome_subunit_indices].min(1).reshape(-1, 1)
	complex_counts[:, :ribosome_subunit_indices] -= inactive_ribosomes  # remove fraction of inactive ribosomes
	rna_counts += complex_counts @ rna_stoich
	protein_counts += complex_counts @ protein_stoich

	# Calculate mass of RNA/protein/ribosomes
	rna_mass = rna_counts @ mw_rnas
	protein_mass = protein_counts @ mw_proteins
	ribosome_mass = (active_ribosome_counts + inactive_ribosomes.squeeze()) * mw_ribosome
	rna_in_ribosome_mass = ribosome_fraction_rna * ribosome_mass
	protein_in_ribosome_mass = ribosome_fraction_protein * ribosome_mass
	total_mass = rna_mass + protein_mass + enzyme_mass + rna_in_ribosome_mass + protein_in_ribosome_mass

	# Excess mass for ribosomal RNA/protein
	excess_rna = rna_mass / ribosome_mass
	excess_protein = protein_mass / ribosome_mass

	# Fractions by synthesis components (ribosome and amino acid synthesis enzymes)
	rna_fraction = (rna_mass + rna_in_ribosome_mass) / total_mass
	protein_fraction = (protein_mass + protein_in_ribosome_mass) / total_mass
	enzyme_fraction = enzyme_mass / total_mass

	# Fractions by protein
	rprotein_protein_fraction = (protein_mass + protein_in_ribosome_mass) / total_protein_mass
	enzyme_protein_fraction = enzyme_mass / total_protein_mass

	excess = np.vstack((excess_rna, excess_protein)).T
	synth_fractions = np.vstack((rna_fraction, protein_fraction, enzyme_fraction)).T
	protein_fractions = np.vstack((rprotein_protein_fraction, enzyme_protein_fraction)).T
	return excess, synth_fractions, protein_fractions


class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)

		# Listeners used
		main_reader = TableReader(os.path.join(simOutDir, 'Main'))

		# Load data
		sim_time = main_reader.readColumn('time')[1:] / 3600
		excess, synth_fractions, protein_fractions = calculate_ribosome_excesses(
			sim_data, [os.path.dirname(simOutDir)])

		# Plot data
		_, (excess_ax, synth_ax, protein_ax) = plt.subplots(3, 1, figsize=(5, 15))

		## Plot excess ribosomes
		excess_ax.plot(sim_time, excess)
		excess_ax.legend(['rRNA', 'rProtein'], fontsize=6, frameon=False)
		excess_ax.set_ylabel('Mass fraction excess of ribosome components')
		self.remove_border(excess_ax)

		## Plot synthesis fractions
		synth_ax.plot(sim_time, synth_fractions)
		synth_ax.legend(['rRNA', 'rProtein', 'Enzymes'], fontsize=6, frameon=False)
		synth_ax.set_ylabel('Mass fractions of synthesis components')
		self.remove_border(synth_ax)

		## Plot protein fractions
		protein_ax.plot(sim_time, protein_fractions)
		protein_ax.legend(['rProtein', 'Enzymes'], fontsize=6, frameon=False)
		protein_ax.set_xlabel('Time (hr)')
		protein_ax.set_ylabel('Mass fractions of all protein')
		self.remove_border(protein_ax)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
