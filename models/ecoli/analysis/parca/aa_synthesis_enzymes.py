"""
Compare expected amino acid synthesis enzyme counts at initialization to
validation data.

TODO:
	- Add subplot for with_aa condition
	- Consider other proteomics datasets
	- Add similar cohort plot for counts seen in sims
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import parcaAnalysisPlot
from reconstruction.ecoli.initialization import create_bulk_container
from wholecell.analysis.analysis_tools import exportFigure


def create_container(sim_data, stoich, complex_ids, subunit_ids, molecule_ids, **options):
	container = create_bulk_container(sim_data, n_seeds=5, **options)
	decomplex_counts = stoich @ container.counts(complex_ids)
	container.countsDec(decomplex_counts, subunit_ids)
	return container.counts(molecule_ids)


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validation_data_file, 'rb') as f:
			validation_data = pickle.load(f)

		metabolism = sim_data.process.metabolism
		complexation = sim_data.process.complexation
		replication = sim_data.process.replication
		translation = sim_data.process.translation

		synthesis_enzymes = metabolism.aa_enzymes
		synthesis_monomers = [
			subunit
			for enzyme in synthesis_enzymes
			for subunit in complexation.get_monomers(enzyme)['subunitIds']
			]

		monomer_to_rna = {monomer['id']: monomer['rna_id'][:-3] for monomer in translation.monomer_data}
		rna_to_symbol = {gene['rna_id']: gene['symbol'] for gene in replication.gene_data}
		synthesis_symbols = [rna_to_symbol[monomer_to_rna[monomer]] for monomer in synthesis_monomers]

		complexation_stoich = complexation.stoich_matrix_monomers()

		wcm_basal_counts_ppgpp_atten = create_container(sim_data, complexation_stoich,
			complexation.ids_complexes, complexation.molecule_names, synthesis_monomers)
		wcm_basal_counts = create_container(sim_data, complexation_stoich,
			complexation.ids_complexes, complexation.molecule_names, synthesis_monomers,
			ppgpp_regulation=False, trna_attenuation=False)
		validation_glucose = dict(zip(validation_data.protein.schmidt2015Data['monomerId'],
			validation_data.protein.schmidt2015Data['glucoseCounts']))
		val_basal_counts = np.array([validation_glucose.get(m, 0) for m in synthesis_monomers])

		plt.figure()

		# Plot data and diagonal reference lines
		max_counts = max(wcm_basal_counts.max(), val_basal_counts.max())
		plt.loglog(val_basal_counts, wcm_basal_counts, 'o', alpha=0.5, label='No ppGpp or attenuation')
		plt.loglog(val_basal_counts, wcm_basal_counts_ppgpp_atten, 'o', alpha=0.5, label='With ppGpp and attenuation')
		plt.loglog([1, max_counts], [1, max_counts], '--k')
		plt.loglog([10, max_counts], [1, max_counts / 10], '--k', linewidth=1)
		plt.loglog([1, max_counts / 10], [10, max_counts], '--k', linewidth=1)

		# Label points with gene symbols
		for x, y, label in zip(val_basal_counts, wcm_basal_counts, synthesis_symbols):
			plt.text(x, 1.2*y, label, ha='center', fontsize=6)

		plt.xlabel('Validation counts')
		plt.ylabel('Model predicted counts')
		self.remove_border()
		plt.legend(fontsize=8, frameon=False)

		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
