"""
Compare the directionality of expression changes for all proteins and groups
related to growth in the model to validation.

TODO:
	plot positive/negative directionality matches
	show absolute number (matched/unmatched) instead of fraction of matches
	scatter plots for each set of regulation options
	compare more environments than just minimal to rich
	shares a lot of code with the variant analysis plot of the same name
"""

import pickle

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

from models.ecoli.analysis import parcaAnalysisPlot
from reconstruction.ecoli.initialization import create_bulk_container
from wholecell.analysis.analysis_tools import exportFigure


def get_mw(mw, molecules):
	return np.array([mw.get(mol) for mol in molecules])

def get_monomers(molecules, get_stoich):
	return [monomer for mol in molecules for monomer in get_stoich(mol)['subunitIds']]

def get_container_fractions(container, molecules, mw):
	total_mass = container.counts(mw.keys()) @ np.array(list(mw.values()))
	return [container.counts(mol) * get_mw(mw, mol) / total_mass for mol in molecules]

def get_validation_fractions(counts, molecules, mw):
	total_mass = np.sum([count * mw.get(mol, 0) for mol, count in counts.items()])
	return [np.array([counts.get(mol, 0) * mw.get(mol, 0) / total_mass for mol in molecule_group]) for molecule_group in molecules]

def compare_fractions(condition1, condition2):
	with np.errstate(divide='ignore', invalid='ignore'):
		fc = [np.log2(c1 / c2) for c1, c2 in zip(condition1, condition2)]
	return fc

def compare_to_validation(parca, validation):
	matches = [(np.sign(p) == np.sign(v))[np.isfinite(p) & np.isfinite(v)] for p, v in zip(parca, validation)]
	return [match.sum() / match.shape[0] for match in matches]


class Plot(parcaAnalysisPlot.ParcaAnalysisPlot):
	def do_plot(self, input_dir, plot_out_dir, plot_out_filename, sim_data_file, validation_data_file, metadata):
		with open(sim_data_file, 'rb') as f:
			sim_data = pickle.load(f)
		with open(validation_data_file, 'rb') as f:
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
		option_choices = {
			'Both': {'ppgpp_regulation': True, 'trna_attenuation': True},
			'Only ppGpp': {'ppgpp_regulation': True, 'trna_attenuation': False},
			'Only attenuation': {'ppgpp_regulation': False, 'trna_attenuation': True},
			'Neither': {'ppgpp_regulation': False, 'trna_attenuation': False},
			}
		parca_compare = {}
		for label, options in option_choices.items():
			rich_container = create_bulk_container(sim_data, condition='with_aa', form_complexes=False, **options)
			basal_container = create_bulk_container(sim_data, form_complexes=False, **options)
			rich_fractions = get_container_fractions(rich_container, monomer_ids, mw)
			basal_fractions = get_container_fractions(basal_container, monomer_ids, mw)
			parca_compare[label] = compare_fractions(rich_fractions, basal_fractions)

		# Validation mass fractions
		rich_fractions_validation = get_validation_fractions(rich_validation, monomer_ids, mw)
		basal_fractions_validation = get_validation_fractions(basal_validation, monomer_ids, mw)
		validation_compare = compare_fractions(rich_fractions_validation, basal_fractions_validation)

		_, (bar_ax, scat_ax) = plt.subplots(2, 1, figsize=(5, 10))

		# Plot bar for fraction of matches between parca and validation changes
		x = np.arange(len(group_labels))
		width = 0.8 / len(parca_compare)
		offset = width / 2 - 0.4
		for i, (label, parca) in enumerate(parca_compare.items()):
			comparison = compare_to_validation(parca, validation_compare)
			bar_ax.bar(x + offset + width * i, comparison, width, label=label)
		bar_ax.set_xticks(x)
		bar_ax.set_xticklabels(group_labels)
		bar_ax.tick_params(labelsize=6)
		bar_ax.set_ylim([0, 1])
		self.remove_border(bar_ax)
		bar_ax.legend(fontsize=6, frameon=False)
		bar_ax.set_ylabel('Fraction direction matches validation', fontsize=8)

		# Plot scatter plot of changes in the parca vs validation
		for val, parca, label in zip(validation_compare, parca_compare['Both'], group_labels):
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
		exportFigure(plt, plot_out_dir, plot_out_filename, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
