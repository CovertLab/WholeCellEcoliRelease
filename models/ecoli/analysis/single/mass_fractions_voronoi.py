"""
Plot the Voronoi diagram of mass fractions

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 09/27/2019
"""

from __future__ import absolute_import, division, print_function
import os
import cPickle
import numpy as np
from matplotlib import pyplot as plt
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import singleAnalysisPlot
from wholecell.utils import units
from wholecell.utils.voronoi_plot_main import VoronoiMaster

SEED = 0 # random seed

class Plot(singleAnalysisPlot.SingleAnalysisPlot):
	def do_plot(self, simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)

		np.random.seed(SEED)

		# Load data
		mass = TableReader(os.path.join(simOutDir, "Mass"))
		unique_molecule_counts = TableReader(os.path.join(simOutDir, "UniqueMoleculeCounts"))
		bulk_molecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulk_molecule_counts = bulk_molecules.readColumn("counts")
		bulk_molecule_idx = {name: i for i, name in enumerate(bulk_molecules.readAttribute("objectNames"))}
		nAvogadro = sim_data.constants.nAvogadro

		# Get the stoichiometry matrix of 50s & 30s subunits
		ribosome_50s_subunits = sim_data.process.complexation.getMonomers(
			sim_data.moleculeIds.s50_fullComplex)
		ribosome_30s_subunits = sim_data.process.complexation.getMonomers(
			sim_data.moleculeIds.s30_fullComplex)
		ribosome_composition_count = np.hstack((ribosome_50s_subunits["subunitStoich"], 
			ribosome_30s_subunits["subunitStoich"]))

		# Count 50s & 30s subunits
		complex_indexes = bulk_molecule_idx[sim_data.moleculeIds.s50_fullComplex]
		complex_50s_counts = bulk_molecule_counts[:, complex_indexes]
		complex_indexes = bulk_molecule_idx[sim_data.moleculeIds.s30_fullComplex]
		complex_30s_counts = bulk_molecule_counts[:, complex_indexes]

		n_50s_subunits = np.outer(complex_50s_counts, ribosome_50s_subunits["subunitStoich"])
		n_30s_subunits = np.outer(complex_30s_counts, ribosome_30s_subunits["subunitStoich"])
		n_ribosome_subunits = np.hstack((n_50s_subunits, n_30s_subunits))

		# Count active ribosomes
		ribosome_subunit_ids = (ribosome_50s_subunits["subunitIds"].tolist() 
			+ ribosome_30s_subunits["subunitIds"].tolist())
		unique_molecule_ids = unique_molecule_counts.readAttribute("objectNames")
		ribosome_idx = unique_molecule_ids.index("active_ribosome")
		n_active_ribosome = unique_molecule_counts.readColumn("uniqueMoleculeCounts")[:, ribosome_idx]
		n_active_ribosome_subunits = np.outer(n_active_ribosome, ribosome_composition_count)

		def find_mass_rRna(rRna_id):
			# Count free rRNAs
			rRna_ids = getattr(sim_data.moleculeGroups, str(rRna_id))
			rRna_indexes = np.array(
				[bulk_molecule_idx[rRna] for rRna in rRna_ids])
			free_rRna_counts = bulk_molecule_counts[:, rRna_indexes]

			# Include the rRNAs already in active ribosomes & ribosome subunits
			total_rRna_counts = free_rRna_counts.astype(float)
			for i, subunit_id in enumerate(ribosome_subunit_ids):
				if subunit_id in rRna_ids:
					rRna_index = rRna_ids.index(subunit_id)
					total_rRna_counts[:, rRna_index] += n_active_ribosome_subunits[:, i]
					total_rRna_counts[:, rRna_index] += n_ribosome_subunits[:, i]

			# Load molecular weights and convert to mass
			rRna_mw = sim_data.getter.getMass(rRna_ids)
			return (units.dot(total_rRna_counts, rRna_mw)/nAvogadro).asNumber(units.fg)

		rRna_16s = find_mass_rRna('s30_16sRRNA')
		rRna_23s = find_mass_rRna('s50_23sRRNA')
		rRna_5s = find_mass_rRna('s50_5sRRNA')

		# lipids and polyamines
		def find_mass_molecule_group(group_id):
			temp_ids = getattr(sim_data.moleculeGroups, str(group_id))
			temp_indexes = np.array([bulk_molecule_idx[temp] for temp in temp_ids])
			temp_counts = bulk_molecule_counts[:, temp_indexes]
			temp_mw = sim_data.getter.getMass(temp_ids)
			return (units.dot(temp_counts, temp_mw)/nAvogadro).asNumber(units.fg)

		lipid = find_mass_molecule_group('lipids')
		polyamines = find_mass_molecule_group('polyamines')

		# LPS, murein, and glycogen
		def find_mass_single_molecule(molecule_id):
			temp_id = getattr(sim_data.moleculeIds, str(molecule_id))
			temp_index = bulk_molecule_idx[temp_id]
			temp_counts = bulk_molecule_counts[:, temp_index]
			temp_mw = sim_data.getter.getMass([temp_id])
			return (units.multiply(temp_counts, temp_mw)/nAvogadro).asNumber(units.fg)

		lps = find_mass_single_molecule('LPS')
		murein = find_mass_single_molecule('murein')
		glycogen = find_mass_single_molecule('glycogen')

		# other cell components
		protein = mass.readColumn("proteinMass")
		rna = mass.readColumn("rnaMass")
		tRna = mass.readColumn("tRnaMass")
		rRna = mass.readColumn("rRnaMass")
		assert np.allclose(rRna, (rRna_16s + rRna_23s + rRna_5s), rtol=1e-10)
		mRna = mass.readColumn("mRnaMass")
		miscRna = rna - (tRna + rRna + mRna)
		dna = mass.readColumn("dnaMass")
		smallMolecules = mass.readColumn("smallMoleculeMass")
		metabolites = smallMolecules - (lipid + lps + murein + polyamines + glycogen)

		# create dictionary
		dic_initial = {
			'nucleic_acid': {
				'DNA': dna[0],
				'mRNA': mRna[0],
				'miscRNA': miscRna[0],
				'rRNA': {
					'16srRNA': rRna_16s[0],
					'23srRNA': rRna_23s[0],
					'5srRNA': rRna_5s[0],
				},
				'tRNA': tRna[0],
			},
			'metabolites': {
				'LPS': lps[0],
				'glycogen': glycogen[0],
				'lipid': lipid[0],
				'metabolites': metabolites[0],
				'peptidoglycan': murein[0],
				'polyamines': polyamines[0],
			},
			'protein': protein[0],
		}
		dic_final = {
			'nucleic_acid': {
				'DNA': dna[-1],
				'mRNA': mRna[-1],
				'miscRNA': miscRna[-1],
				'rRNA': {
					'16srRNA': rRna_16s[-1],
					'23srRNA': rRna_23s[-1],
					'5srRNA': rRna_5s[-1],
				},
				'tRNA': tRna[-1],
			},
			'metabolites': {
				'LPS': lps[-1],
				'glycogen': glycogen[-1],
				'lipid': lipid[-1],
				'metabolites': metabolites[-1],
				'peptidoglycan': murein[-1],
				'polyamines': polyamines[-1],
			},
			'protein': protein[-1],
		}

		# create the plot
		vm = VoronoiMaster()
		vm.plot([[dic_initial, dic_final]],
				title = [["Initial biomass components", "Final biomass components"]],
				ax_shape = (1, 2), chained = True)

		# save the plot
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close("all")
		
if __name__ == "__main__":
	Plot().cli()
