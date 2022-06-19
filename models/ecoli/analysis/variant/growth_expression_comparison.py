"""
Compare the directionality of expression changes for all proteins and groups
related to growth in the model to validation.

TODO:
	shares a lot of code with the parca analysis plot of the same name
"""

import pickle
import os

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import exportFigure, read_stacked_columns
from wholecell.io.tablereader import TableReader


def read_counts(cell_paths):
	ids = TableReader(os.path.join(cell_paths[0], 'simOut', 'MonomerCounts')).readAttribute('monomerIds')
	counts = read_stacked_columns(cell_paths, 'MonomerCounts', 'monomerCounts').T
	return dict(zip(ids, counts))

def get_mw(mw, molecules):
	return np.array([mw.get(mol) for mol in molecules])

def get_monomers(molecules, get_stoich):
	return [monomer for mol in molecules for monomer in get_stoich(mol)['subunitIds']]

def get_sim_fractions(counts, molecules, mw):
	total_mass = np.sum([count * mw.get(mol, 0) for mol, count in counts.items()], axis=0)
	return [
		np.array([counts.get(mol, np.zeros_like(total_mass)) * mw.get(mol, 0) / total_mass for mol in molecule_group]).mean(1)
		for molecule_group in molecules
		]

def get_validation_fractions(counts, molecules, mw):
	total_mass = np.sum([count * mw.get(mol, 0) for mol, count in counts.items()])
	return [
		np.array([counts.get(mol, 0) * mw.get(mol, 0) / total_mass for mol in molecule_group])
		for molecule_group in molecules
		]

def compare_fractions(condition1, condition2):
	with np.errstate(divide='ignore', invalid='ignore'):
		fc = [np.log2(c1 / c2) for c1, c2 in zip(condition1, condition2)]
	return fc

def compare_to_validation(parca, validation):
	matches = [(np.sign(p) == np.sign(v))[np.isfinite(p) & np.isfinite(v)] for p, v in zip(parca, validation)]
	return [match.sum() / match.shape[0] for match in matches]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		variants = self.ap.get_variants()

		minimal_sim = None
		rich_sim = None
		for variant in variants:
			with open(self.ap.get_variant_kb(variant), 'rb') as f:
				variant_sim_data = pickle.load(f)

			cell_paths = self.ap.get_cells(variant=[variant])
			if variant_sim_data.condition == 'basal':
				minimal_sim = read_counts(cell_paths)
			elif variant_sim_data.condition == 'with_aa':
				rich_sim = read_counts(cell_paths)

		if minimal_sim is None or rich_sim is None:
			print(f'Do not have a minimal and rich variant for {plotOutFileName} - skipping...')
			return

		with open(simDataFile, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = pickle.load(f)

		# sim_data attributes used
		metabolism = sim_data.process.metabolism
		transcription = sim_data.process.transcription
		delta_prob = sim_data.process.transcription_regulation.delta_prob
		translation = sim_data.process.translation
		mol_ids = sim_data.molecule_ids
		get_stoich = sim_data.process.complexation.get_monomers
		aa_enzymes = [e for e in metabolism.aa_enzymes if e != 'PUTA-CPLXBND[c]']  # PutA is already accounted for as PUTA-CPLX and no easy way to get monomers from both complexes and equilibrium molecules (like PUTA-CPLXBND)

		# Get regulated subgroup
		tf_regulated_indices = np.unique([
			i for i, v in zip(delta_prob['deltaI'], delta_prob['deltaV'])
			if v != 0
			])
		cistron_id_to_idx = {
			cistron: i
			for i, cistron in enumerate(transcription.cistron_data['id'])
			}
		ppgpp_regulated_indices = np.array([
			cistron_id_to_idx[cistron]
			for cistron in transcription.ppgpp_regulated_genes
			])
		regulated_cistrons = np.zeros(len(transcription.cistron_data), bool)
		regulated_cistrons[tf_regulated_indices] = True
		regulated_cistrons[ppgpp_regulated_indices] = True
		regulated_cistrons[transcription.attenuated_rna_indices] = True
		regulated_monomers = regulated_cistrons[sim_data.relation.cistron_to_monomer_mapping]

		# Validation data
		rich_validation = {}
		basal_validation = {}
		for p in validation_data.protein.li_2014:
			monomer = p['monomer']
			rich_validation[monomer] = p['rich_rate']
			basal_validation[monomer] = p['minimal_rate']

		mw = {monomer['id']: monomer['mw'] for monomer in translation.monomer_data}

		# Select molecule groups of interest
		monomer_ids = (
			translation.monomer_data['id'],
			translation.monomer_data['id'][regulated_monomers],
			get_monomers(aa_enzymes, get_stoich),
			get_monomers(transcription.synthetase_names, get_stoich),
			get_monomers([mol_ids.RelA, mol_ids.SpoT], get_stoich),
		)
		group_labels = [
			'All',
			'Regulated',
			'AA enzymes',
			'Synthetases',
			'ppGpp molecules',
		]

		# Expected bulk containers in different conditions
		rich_fractions = get_sim_fractions(rich_sim, monomer_ids, mw)
		basal_fractions = get_sim_fractions(minimal_sim, monomer_ids, mw)
		parca_compare = compare_fractions(rich_fractions, basal_fractions)

		# Validation mass fractions
		rich_fractions_validation = get_validation_fractions(rich_validation, monomer_ids, mw)
		basal_fractions_validation = get_validation_fractions(basal_validation, monomer_ids, mw)
		validation_compare = compare_fractions(rich_fractions_validation, basal_fractions_validation)

		comparison = compare_to_validation(parca_compare, validation_compare)

		_, (bar_ax, scat_ax) = plt.subplots(2, 1, figsize=(5, 10))

		# Plot bar for fraction of matches between sim and validation changes
		bar_ax.bar(group_labels, comparison)
		bar_ax.tick_params(labelsize=6)
		bar_ax.set_ylim([0, 1])
		self.remove_border(bar_ax)
		bar_ax.set_ylabel('Fraction direction matches validation', fontsize=8)

		# Plot scatter plot of changes in the sim vs validation
		for val, parca, label in zip(validation_compare, parca_compare, group_labels):
			mask = np.isfinite(val) & np.isfinite(parca)
			n = np.sum(mask)
			_, r = stats.pearsonr(val[mask], parca[mask])
			if n > 200:
				options = {'alpha': 0.2, 'markersize': 1}
			else:
				options = {'alpha': 0.5, 'markersize': 2}
			scat_ax.plot(val, parca, 'o', label=f'{label} (r={r:.3f}, n={n})', **options)
		scat_ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
		scat_ax.axvline(0, color='k', linestyle='--', linewidth=0.5)
		scat_ax.set_xlabel('Validation log2 fold change\nfrom minimal to rich', fontsize=8)
		scat_ax.set_ylabel('Sim log2 fold change\nfrom minimal to rich', fontsize=8)
		scat_ax.tick_params(labelsize=6)
		self.remove_border(scat_ax)
		scat_ax.legend(fontsize=6, frameon=False)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
