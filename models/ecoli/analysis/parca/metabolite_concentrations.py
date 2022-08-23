"""
Compare metabolite concentrations from different datasets.

TODO:
	include species/conc from getBiomassAsConcentrations
"""

from __future__ import absolute_import, division, print_function

import os
from typing import cast, Dict, List

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy import stats
from six.moves import cPickle, range

from models.ecoli.analysis import parcaAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import constants, units
from six.moves import zip


CONC_UNITS = units.mmol / units.L
INDEX_KEY = 'index'
CONC_KEY = 'conc'
WCM_KEY = 'WCM'
SANDER_KEY = 'Sander Concentration'


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(os.path.join(input_dir, constants.SERIALIZED_RAW_DATA), 'rb') as f:
			raw_data = cPickle.load(f)
		with open(sim_data_file, 'rb') as f:
			sim_data = cPickle.load(f)

		# Extract raw concentrations
		concentrations = {}  # type: Dict[str, Dict[str, List]]
		metabolites = []
		index = 0
		for row in raw_data.metabolite_concentrations:
			metabolites.append(row['Metabolite'])
			for source, conc in row.items():
				if source == 'Metabolite':
					continue

				conc = conc.asNumber(CONC_UNITS)
				if np.isfinite(conc):
					if source not in concentrations:
						concentrations[source] = {
							INDEX_KEY: [],
							CONC_KEY: [],
							}
					concentrations[source][INDEX_KEY].append(index)
					concentrations[source][CONC_KEY].append(conc)
			index += 1

		# Extract sim_data concentrations
		model_conc = np.array([
			sim_data.process.metabolism.conc_dict[met + '[c]'].asNumber(CONC_UNITS)
			for met in metabolites
			])

		# Sort to determine x position
		sorted_idx = np.argsort(model_conc)[::-1]
		sorted_mapping = {idx: i for i, idx in enumerate(sorted_idx)}
		x = list(range(len(metabolites)))

		# Create plot
		plt.figure(figsize=(15, 10))

		## Plot raw sources
		for source, data in concentrations.items():
			new_idx = [sorted_mapping[i] for i in data[INDEX_KEY]]
			plt.semilogy(new_idx, data[CONC_KEY], 'o', alpha=0.5, markersize=4,
				markeredgewidth=0, label=source)

		## Plot wholecell concentration
		plt.semilogy(x, model_conc[sorted_idx],
			'_k', alpha=0.8, label='WCM')

		## Draw reference lines
		for pos in x[::5]:
			plt.axvline(pos - 0.5, color='k', alpha=0.2, linewidth=0.5)

		## Formatting
		plt.xticks(x, np.array(metabolites)[sorted_idx], rotation=45, ha='right', size=6)
		plt.ylabel('Concentration (mM)')
		plt.legend()

		## Save figure
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')

		# Compare sources against each other
		plt.figure(figsize=(10, 10))
		concentrations[WCM_KEY] = {}
		concentrations[WCM_KEY][INDEX_KEY] = x
		concentrations[WCM_KEY][CONC_KEY] = cast(list, model_conc.tolist())
		sources = sorted(concentrations.keys())
		n_sources = len(sources)

		## Subplot for each comparison
		gs = gridspec.GridSpec(n_sources, n_sources)
		for i, source1 in enumerate(sources):
			for j, source2 in enumerate(sources):
				# Extract common concentrations
				idx1 = set(concentrations[source1][INDEX_KEY])
				idx2 = set(concentrations[source2][INDEX_KEY])
				conc1 = [c for c, idx in zip(concentrations[source1][CONC_KEY], concentrations[source1][INDEX_KEY]) if idx in idx2]
				conc2 = [c for c, idx in zip(concentrations[source2][CONC_KEY], concentrations[source2][INDEX_KEY]) if idx in idx1]
				conc_range = [np.min([conc1, conc2]), np.max([conc1, conc2])]
				r, p = stats.pearsonr(np.log(conc1), np.log(conc2))

				# Plot data
				plt.subplot(gs[j, i])
				plt.loglog(conc1, conc2, 'o', alpha=0.5)
				plt.loglog(conc_range, conc_range, 'k--')
				plt.title(r'$R^2$={:.2f}, p={:.1g}'.format(r, p), fontsize=8)

				# Only show axis labels on edge
				if j == n_sources - 1:
					plt.xlabel(source1 + '\n(mM)', fontsize=8)
				if i == 0:
					plt.ylabel(source2 + '\n(mM)', fontsize=8)

		## Save figure
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename + '_sources', metadata)
		plt.close('all')

		# Plot amino acids in Sander paper
		metabolite_idx = {m: i for i, m in enumerate(metabolites)}
		allosteric_aas = ['ARG', 'TRP', 'HIS', 'LEU', 'ILE', 'THR', 'PRO']
		allosteric_idx = [metabolite_idx[aa] for aa in allosteric_aas]

		plt.figure()
		for source, data in concentrations.items():
			conc = {i: c for i, c in zip(data[INDEX_KEY], data[CONC_KEY])}
			x = []
			y = []
			for i, idx in enumerate(allosteric_idx):
				if idx in conc:
					x.append(i)
					y.append(conc[idx])

			marker = 'o'
			options = dict(alpha=0.5, markersize=6, markeredgewidth=0)
			if source == WCM_KEY:
				marker = '_k'
				options = dict(alpha=0.8)
			elif source == SANDER_KEY:
				marker = 'X'
				options['markersize'] = 8

			plt.semilogy(x, y, marker, label=source, **options)

		plt.xticks(range(len(allosteric_aas)), allosteric_aas)
		plt.ylabel('Concentration (mM)')
		plt.legend()

		## Save figure
		plt.tight_layout()
		exportFigure(plt, plot_out_dir, plot_out_filename + '_sander', metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
