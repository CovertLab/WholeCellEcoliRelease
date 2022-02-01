"""
Compare amino acid concentrations across variants similar to data from
Sander et al. Allosteric Feedback Inhibition Enables Robust Amino Acid
Biosynthesis in E. coli by Enforcing Enzyme Overabundance. 2019. Fig 1B.

Associated variant to run:
	remove_aa_inhibition

TODO:
	- add reference conc from Sander et al
"""

import os

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.sim.variants.remove_aa_inhibition import AA_TO_ENZYME, get_aa_and_ki_factor
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


WILDTYPE = 'wt'
HEATMAP_ROWS = ['hisG', 'proB', 'trpE', 'argA', WILDTYPE, 'ilvA', 'leuA', 'thrA']
HEATMAP_COLS = [
	('gln', ['GLN[c]']),
	('glu', ['GLT[c]']),
	('arg', ['ARG[c]']),
	('pro', ['PRO[c]']),
	('asp', ['L-ASPARTATE[c]']),
	('asn', ['ASN[c]']),
	('lys', ['LYS[c]']),
	('met', ['MET[c]']),
	('thr', ['THR[c]']),
	('ile/leu', ['ILE[c]', 'LEU[c]']),
	('val', ['VAL[c]']),
	('ala', ['L-ALPHA-ALANINE[c]']),
	('ser', ['SER[c]']),
	('gly', ['GLY[c]']),
	('his', ['HIS[c]']),
	('phe', ['PHE[c]']),
	('trp', ['TRP[c]']),
	('tyr', ['TYR[c]']),
	('ile', ['ILE[c]']),
	('leu', ['LEU[c]']),
]


def subplot(gs, aa_idx, control, control_variance, data, variance, xlabels, title, amino_acids, factors=None, legend=None):
	ax = plt.subplot(gs)
	_, n_variants, n_factors = data.shape
	x = np.array(list(range(n_variants))) + 1
	width = 0.8 / n_factors
	offsets = np.array(list(range(n_factors))) * width - 0.4 + width/2
	idx = np.array([aa_idx[aa] for aa in amino_acids])

	# Plot an amino acid concentration (or sum of multiple amino acids)
	# for each mutant (variant)
	error_kw = dict(elinewidth=2*width)
	ax.bar(0, control[idx].sum(), width, label=WILDTYPE,
		yerr=np.sqrt(control_variance[idx].sum()), error_kw=error_kw)
	for i, offset in enumerate(offsets):
		label = None if factors is None else f'KI x {factors[i]}'
		ax.bar(x + offset, data[idx, :, i].sum(axis=0), width, label=label,
			yerr=np.sqrt(variance[idx, :, i].sum(axis=0)), error_kw=error_kw)

	# Format subplot
	ax.set_yscale('log')
	ax.set_xticks([0] + list(x))
	ax.set_xticklabels(xlabels, fontsize=8, rotation=45)
	ax.set_ylabel('Conc (mM)', fontsize=8)
	ax.set_title(title, fontsize=8)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	if legend:
		handles, labels = ax.get_legend_handles_labels()
		legend_ax = plt.subplot(legend, frameon=False)
		legend_ax.axis('off')
		legend_ax.legend(handles, labels, loc='center', frameon=False)

def heatmap(gs, aa_idx, data, enzyme_order):
	ax = plt.subplot(gs)
	n_rows = len(HEATMAP_ROWS)
	n_cols = len(HEATMAP_COLS)

	# Extract concentrations from data
	conc = np.zeros((n_rows, n_cols))
	for row_idx, enzyme in enumerate(HEATMAP_ROWS):
		if enzyme in enzyme_order:
			variant_idx = enzyme_order.index(enzyme)
			for col_idx, cols in enumerate(HEATMAP_COLS):
				for aa in cols[1]:
					conc[row_idx, col_idx] += data[variant_idx, aa_idx[aa]]

	# Scale concentrations for heatmap
	wt_idx = HEATMAP_ROWS.index(WILDTYPE)
	conc /= conc[wt_idx, :]
	conc_fc = np.log2(conc)

	# Plot heatmap
	im = plt.imshow(conc_fc, cmap='RdBu')

	# Show colorbar
	lim = max(-conc_fc.min(), conc_fc.max())
	plt.clim(-lim, lim)
	plt.colorbar(im)

	# Format axes
	plt.xticks(range(n_cols), [col[0] for col in HEATMAP_COLS], rotation=90)
	plt.yticks(range(n_rows), HEATMAP_ROWS)
	ax.set_title('Relative Concentration Change (log2)')


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	_suppress_numpy_warnings = True

	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		variants = self.ap.get_variants()

		aa_ids = sorted({aa for aas in HEATMAP_COLS for aa in aas[1]})
		aa_idx = {aa: i for i, aa in enumerate(aa_ids)}
		n_aas = len(aa_ids)

		control_conc = np.zeros(n_aas)
		control_var = np.zeros(n_aas)
		aa_conc = {}
		aa_var = {}
		aa_variants = []
		ki_factors = []
		for i, variant in enumerate(variants):
			variant_conc = []
			for sim_dir in self.ap.get_cells(variant=[variant], only_successful=True):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				kinetics_reader = TableReader(os.path.join(simOutDir, 'EnzymeKinetics'))

				# Read data
				(aa_counts,) = read_bulk_molecule_counts(simOutDir, aa_ids)
				counts_to_molar = kinetics_reader.readColumn('countsToMolar').reshape(-1, 1)

				# Calculate amino acid concentration at each time step
				variant_conc.append(aa_counts[1:, :]*counts_to_molar[1:])

			# Average concentration over all time steps
			conc = np.vstack(variant_conc)
			conc_mean = conc.mean(axis=0)
			conc_var = conc.std(axis=0)**2
			if variant == 0:
				control_conc = conc_mean
				control_var = conc_var
			else:
				key = get_aa_and_ki_factor(variant)
				aa_variants.append(key[0])
				ki_factors.append(key[1])
				aa_conc[key] = conc_mean
				aa_var[key] = conc_var

		# Process data into a new matrix for easier handling on bar charts
		aa_variants = [aa for aa in AA_TO_ENZYME if aa in aa_variants]
		ki_factors = np.unique(ki_factors)[::-1]
		conc = np.zeros((n_aas, len(aa_variants), len(ki_factors)))
		variance = np.zeros((n_aas, len(aa_variants), len(ki_factors)))
		for i, aa, in enumerate(aa_variants):
			for j, factor in enumerate(ki_factors):
				key = (aa, factor)
				if key in aa_conc:
					conc[:, i, j] = aa_conc[key]
					variance[:, i, j] = aa_var[key]

		# Skip plotting if data does not exist
		if len(aa_variants) == 0 or len(ki_factors) == 0:
			print(f'Not enough data to plot for {plotOutFileName}')
			return

		# xtick labels
		labels = [WILDTYPE] + list(AA_TO_ENZYME.values())

		# Create figure
		plt.figure(figsize=(10, 25))
		gs = gridspec.GridSpec(nrows=8, ncols=3)

		## Plot heatmap for all amino acids only if data exists for all mutants
		if conc.shape[1] == len(AA_TO_ENZYME):
			heatmap_data = np.vstack((control_conc, conc[:, :, 0].T))
			heatmap(gs[:2, :], aa_idx, heatmap_data, labels)
		else:
			print(f'Not enough data to plot for the heatmap in {plotOutFileName}')

		## Plot subplots for each amino acid with only the inhibition removed to match paper
		no_inhibition_conc = conc[:, :, :1]
		no_inhibition_var = variance[:, :, :1]
		subplot(gs[2, 0], aa_idx, control_conc, control_var, no_inhibition_conc,
			no_inhibition_var, labels, 'Arginine', ['ARG[c]'])
		subplot(gs[2, 1], aa_idx, control_conc, control_var, no_inhibition_conc,
			no_inhibition_var, labels, 'Tryptophan', ['TRP[c]'])
		subplot(gs[2, 2], aa_idx, control_conc, control_var, no_inhibition_conc,
			no_inhibition_var, labels, 'Histidine', ['HIS[c]'])
		subplot(gs[3, 0], aa_idx, control_conc, control_var, no_inhibition_conc,
			no_inhibition_var, labels, '(Iso-)leucine', ['ILE[c]', 'LEU[c]'])  # group isoforms like paper
		subplot(gs[3, 1], aa_idx, control_conc, control_var, no_inhibition_conc,
			no_inhibition_var, labels, 'Threonine', ['THR[c]'])
		subplot(gs[3, 2], aa_idx, control_conc, control_var, no_inhibition_conc,
			no_inhibition_var, labels, 'Proline', ['PRO[c]'])
		subplot(gs[4, 0], aa_idx, control_conc, control_var, no_inhibition_conc,
			no_inhibition_var, labels, 'Isoleucine', ['ILE[c]'])  # also plot single isoforms since model allows granularity
		subplot(gs[4, 1], aa_idx, control_conc, control_var, no_inhibition_conc,
			no_inhibition_var, labels, 'Leucine', ['LEU[c]'], factors=ki_factors, legend=gs[4, 2])  # also plot single isoforms since model allows granularity

		## Plot subplots for each amino acid and show results each KI factor
		subplot(gs[5, 0], aa_idx, control_conc, control_var, conc, variance, labels,
			'Arginine', ['ARG[c]'])
		subplot(gs[5, 1], aa_idx, control_conc, control_var, conc, variance, labels,
			'Tryptophan', ['TRP[c]'])
		subplot(gs[5, 2], aa_idx, control_conc, control_var, conc, variance, labels,
			'Histidine', ['HIS[c]'])
		subplot(gs[6, 0], aa_idx, control_conc, control_var, conc, variance, labels,
			'(Iso-)leucine', ['ILE[c]', 'LEU[c]'])  # group isoforms like paper
		subplot(gs[6, 1], aa_idx, control_conc, control_var, conc, variance, labels,
			'Threonine', ['THR[c]'])
		subplot(gs[6, 2], aa_idx, control_conc, control_var, conc, variance, labels,
			'Proline', ['PRO[c]'])
		subplot(gs[7, 0], aa_idx, control_conc, control_var, conc, variance, labels,
			'Isoleucine', ['ILE[c]'])  # also plot single isoforms since model allows granularity
		subplot(gs[7, 1], aa_idx, control_conc, control_var, conc, variance, labels,
			'Leucine', ['LEU[c]'], factors=ki_factors, legend=gs[7, 2])  # also plot single isoforms since model allows granularity

		## Format and save figure
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
