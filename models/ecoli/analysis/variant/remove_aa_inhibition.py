"""
Compare amino acid concentrations across variants similar to data from
Sander et al. Allosteric Feedback Inhibition Enables Robust Amino Acid
Biosynthesis in E. coli by Enforcing Enzyme Overabundance. 2019. Fig 1B.

Associated variant to run:
	remove_aa_inhibition
"""

import os

from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.remove_aa_inhibition import AA_TO_ENZYME
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


def subplot(gs, aa_idx, data, variance, xlabels, title, amino_acids):
	ax = plt.subplot(gs)
	n_variants = data.shape[0]
	x = list(range(n_variants))
	idx = np.array([aa_idx[aa] for aa in amino_acids])

	# Plot an amino acid concentration (or sum of multiple amino acids)
	# for each mutant (variant)
	ax.bar(x, data[:, idx].sum(axis=1), yerr=np.sqrt(variance[:, idx].sum(axis=1)))

	# Format subplot
	if len(xlabels) == n_variants:
		ax.set_xticks(x)
		ax.set_xticklabels(xlabels, fontsize=8, rotation=45)
	ax.set_ylabel('Conc (mM)', fontsize=8)
	ax.set_title(title, fontsize=8)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

def heatmap(gs, aa_idx, data, enzyme_order):
	ax = plt.subplot(gs)
	n_rows = len(HEATMAP_ROWS)
	n_cols = len(HEATMAP_COLS)

	# Extract concentrations from data
	conc = np.zeros((n_rows, n_cols))
	for row_idx, enzyme in enumerate(HEATMAP_ROWS):
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
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		aa_ids = list({aa for aas in HEATMAP_COLS for aa in aas[1]})
		aa_idx = {aa: i for i, aa in enumerate(aa_ids)}

		aa_conc = np.zeros((len(variants), len(aa_ids)))
		aa_var = np.zeros((len(variants), len(aa_ids)))
		for i, variant in enumerate(variants):
			variant_conc = []
			for sim_dir in ap.get_cells(variant=[variant]):
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
			aa_conc[i, :] = conc.mean(axis=0)
			aa_var[i, :] = conc.std(axis=0)**2

		# xtick labels
		labels = [WILDTYPE] + list(AA_TO_ENZYME.values())

		# Create figure
		plt.figure(figsize=(10, 15))
		gs = gridspec.GridSpec(nrows=5, ncols=3)

		## Plot heatmap for all amino acids
		heatmap(gs[:2, :], aa_idx, aa_conc, labels)

		## Plot subplots for each amino acid
		subplot(gs[2, 0], aa_idx, aa_conc, aa_var, labels, 'Arginine', ['ARG[c]'])
		subplot(gs[2, 1], aa_idx, aa_conc, aa_var, labels, 'Tryptophan', ['TRP[c]'])
		subplot(gs[2, 2], aa_idx, aa_conc, aa_var, labels, 'Histidine', ['HIS[c]'])
		subplot(gs[3, 0], aa_idx, aa_conc, aa_var, labels, '(Iso-)leucine', ['ILE[c]', 'LEU[c]'])  # group isoforms like paper
		subplot(gs[3, 1], aa_idx, aa_conc, aa_var, labels, 'Threonine', ['THR[c]'])
		subplot(gs[3, 2], aa_idx, aa_conc, aa_var, labels, 'Proline', ['PRO[c]'])
		subplot(gs[4, 0], aa_idx, aa_conc, aa_var, labels, 'Isoleucine', ['ILE[c]'])  # also plot single isoforms since model allows granularity
		subplot(gs[4, 1], aa_idx, aa_conc, aa_var, labels, 'Leucine', ['LEU[c]'])  # also plot single isoforms since model allows granularity

		## Format and save figure
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
