"""
Compare proteome for proteins with measured degradation rates.
Only runs with two variants (old degradation rates and new degradation rates).

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/30/19
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, filepath


PROTEINS_WITH_HALF_LIFE = [
	'DCUR-MONOMER[c]',
	'DETHIOBIOTIN-SYN-MONOMER[c]',
	'EG10863-MONOMER[c]',
	'CARBPSYN-SMALL[c]',
	'EG10743-MONOMER[c]',
	'GLUTCYSLIG-MONOMER[c]',
	'CDPDIGLYSYN-MONOMER[i]',
	]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		with open(os.path.join(inputDir, 'kb', constants.SERIALIZED_FIT1_FILENAME), 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		expected_n_variants = 2
		n_variants = len(variants)

		if n_variants < expected_n_variants:
			print('This plot only runs for {} variants.'.format(expected_n_variants))
			return

		# IDs for appropriate proteins
		ids_complexation = sim_data.process.complexation.moleculeNames
		ids_complexation_complexes = sim_data.process.complexation.ids_complexes
		ids_equilibrium = sim_data.process.equilibrium.moleculeNames
		ids_equilibrium_complexes = sim_data.process.equilibrium.ids_complexes
		ids_translation = sim_data.process.translation.monomerData['id'].tolist()
		ids_protein = sorted(set(ids_complexation + ids_equilibrium + ids_translation))

		# Stoichiometry matrices
		equil_stoich = sim_data.process.equilibrium.stoichMatrixMonomers()
		complex_stoich = sim_data.process.complexation.stoichMatrixMonomers()

		# Protein container views
		protein_container = BulkObjectsContainer(ids_protein, dtype=np.float64)
		view_complexation = protein_container.countsView(ids_complexation)
		view_complexation_complexes = protein_container.countsView(ids_complexation_complexes)
		view_equilibrium = protein_container.countsView(ids_equilibrium)
		view_equilibrium_complexes = protein_container.countsView(ids_equilibrium_complexes)

		# Load model data
		model_counts = np.zeros((len(PROTEINS_WITH_HALF_LIFE), expected_n_variants))
		model_std = np.zeros((len(PROTEINS_WITH_HALF_LIFE), expected_n_variants))
		for i, variant in enumerate(variants):
			if i >= expected_n_variants:
				print('Skipping variant {} - only runs for {} variants.'.format(variant, expected_n_variants))
				continue

			variant_counts = []
			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, 'simOut')

				# Listeners used
				unique_counts_reader = TableReader(os.path.join(simOutDir, 'UniqueMoleculeCounts'))

				# Account for bulk molecules
				(bulk_counts,) = read_bulk_molecule_counts(simOutDir, ids_protein)
				protein_container.countsIs(bulk_counts.mean(axis=0))

				# Account for unique molecules
				ribosome_index = unique_counts_reader.readAttribute('uniqueMoleculeIds').index('activeRibosome')
				rnap_index = unique_counts_reader.readAttribute('uniqueMoleculeIds').index('activeRnaPoly')
				n_ribosomes = unique_counts_reader.readColumn('uniqueMoleculeCounts')[:, ribosome_index]
				n_rnap = unique_counts_reader.readColumn('uniqueMoleculeCounts')[:, rnap_index]
				protein_container.countsInc(n_ribosomes.mean(), [sim_data.moleculeIds.s30_fullComplex, sim_data.moleculeIds.s50_fullComplex])
				protein_container.countsInc(n_rnap.mean(), [sim_data.moleculeIds.rnapFull])

				# Account for small-molecule bound complexes
				view_equilibrium.countsDec(equil_stoich.dot(view_equilibrium_complexes.counts()))

				# Account for monomers in complexed form
				view_complexation.countsDec(complex_stoich.dot(view_complexation_complexes.counts()))

				variant_counts.append(protein_container.countsView(PROTEINS_WITH_HALF_LIFE).counts())
			model_counts[:, i] = np.mean(variant_counts, axis=0)
			model_std[:, i] = np.std(variant_counts, axis=0)

		# Validation data
		schmidt_ids = {m: i for i, m in enumerate(validation_data.protein.schmidt2015Data['monomerId'])}
		schmidt_counts = validation_data.protein.schmidt2015Data['glucoseCounts']
		validation_counts = np.array([schmidt_counts[schmidt_ids[p]] for p in PROTEINS_WITH_HALF_LIFE])

		# Process data
		model_log_counts = np.log10(model_counts)
		model_log_lower_std = model_log_counts - np.log10(model_counts - model_std)
		model_log_upper_std = np.log10(model_counts + model_std) - model_log_counts
		validation_log_counts = np.log10(validation_counts)
		r_before = stats.pearsonr(validation_log_counts, model_log_counts[:, 0])
		r_after = stats.pearsonr(validation_log_counts, model_log_counts[:, 1])

		# Scatter plot of model vs validation counts
		max_counts = np.ceil(max(validation_log_counts.max(), model_log_upper_std.max()))
		limits = [0, max_counts]
		plt.figure()
		colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

		## Plot data
		for i in range(expected_n_variants):
			plt.errorbar(validation_log_counts, model_log_counts[:, i],
				yerr=np.vstack((model_log_lower_std[:, i], model_log_upper_std[:, i])),
				fmt='o', color=colors[i], ecolor='k', capsize=3, alpha=0.5)
		plt.plot(limits, limits, 'k--', linewidth=0.5, label='_nolegend_')

		## Format axes
		plt.xlabel('Validation Counts\n(log10(counts))')
		plt.ylabel('Average Simulation Counts\n(log10(counts))')
		ax = plt.gca()
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['left'].set_position(('outward', 10))
		ax.spines['bottom'].set_position(('outward', 10))
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		ax.yaxis.set_major_locator(MaxNLocator(integer=True))

		## Add legend
		legend_text = [
			'Before: r={:.2f}, p={:.3f}'.format(r_before[0], r_before[1]),
			'After: r={:.2f}, p={:.3f}'.format(r_after[0], r_after[1]),
			]
		plt.legend(legend_text, frameon=False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
