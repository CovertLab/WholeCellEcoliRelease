"""
Plot the Voronoi diagram of mass fractions

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 09/27/2019
"""

from __future__ import absolute_import, division, print_function
import os
from six.moves import cPickle
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
		bulk_molecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
		bulk_molecule_counts = bulk_molecules.readColumn("counts")
		bulk_molecule_idx = {name: i for i, name in enumerate(bulk_molecules.readAttribute("objectNames"))}
		nAvogadro = sim_data.constants.n_Avogadro

		# lipids and polyamines
		def find_mass_molecule_group(group_id):
			temp_ids = getattr(sim_data.molecule_groups, str(group_id))
			temp_indexes = np.array([bulk_molecule_idx[temp] for temp in temp_ids])
			temp_counts = bulk_molecule_counts[:, temp_indexes]
			temp_mw = sim_data.getter.get_mass(temp_ids)
			return (units.dot(temp_counts, temp_mw)/nAvogadro).asNumber(units.fg)

		lipid = find_mass_molecule_group('lipids')
		polyamines = find_mass_molecule_group('polyamines')

		# LPS, murein, and glycogen
		def find_mass_single_molecule(molecule_id):
			temp_id = getattr(sim_data.molecule_ids, str(molecule_id))
			temp_index = bulk_molecule_idx[temp_id]
			temp_counts = bulk_molecule_counts[:, temp_index]
			temp_mw = sim_data.getter.get_mass([temp_id])
			return (units.multiply(temp_counts, temp_mw)/nAvogadro).asNumber(units.fg)

		lps = find_mass_single_molecule('LPS')
		murein = find_mass_single_molecule('murein')
		glycogen = find_mass_single_molecule('glycogen')

		# other cell components
		protein = mass.readColumn("proteinMass")
		rna = mass.readColumn("rnaMass")
		tRna = mass.readColumn("tRnaMass")
		rRna = mass.readColumn("rRnaMass")
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
				'rRNA': rRna[0],
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
				'rRNA': rRna[-1],
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
