"""
Analysis for time_step variant to show impact of a maximum time step limit
for processes.
"""

from __future__ import absolute_import, division, print_function

from six.moves import cPickle
import os

from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.time_step import TIME_STEP_FACTOR
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader


def remove_border(ax):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.set_xticks([])
	ax.tick_params(axis='y', labelsize=6)

def show_xaxis(ax, x_vals):
	ax.spines['bottom'].set_visible(True)
	ax.set_xticks(x_vals)
	ax.set_xticklabels(np.array(TIME_STEP_FACTOR)[x_vals], rotation=45, fontsize=6)
	ax.set_xlabel('Fraction of max time step', fontsize=8)

def plot_bar(gs, x, y, y_label, show_x=False):
	ax = plt.subplot(gs)
	ax.bar(x, y)
	ax.set_ylabel(y_label, fontsize=8)
	remove_border(ax)
	if show_x:
		show_xaxis(ax, x)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)

		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		inactive_rnap_id = [sim_data.molecule_ids.full_RNAP]
		ribosome_subunit_ids = [
			sim_data.molecule_ids.s50_full_complex,
			sim_data.molecule_ids.s30_full_complex,
			]

		all_time_steps = np.zeros(n_variants)
		all_doubling_times = np.zeros(n_variants)
		all_execution_times = np.zeros(n_variants)
		all_rnap_elong_rates = np.zeros(n_variants)
		all_rnap_fraction_active = np.zeros(n_variants)
		all_rib_elong_rates = np.zeros(n_variants)
		all_rib_fraction_active = np.zeros(n_variants)
		all_protein_fold_changes = np.zeros(n_variants)
		all_dna_fold_changes = np.zeros(n_variants)
		all_sm_fold_changes = np.zeros(n_variants)
		all_trna_fold_changes = np.zeros(n_variants)
		all_rrna_fold_changes = np.zeros(n_variants)
		all_mrna_fold_changes = np.zeros(n_variants)
		x_vals = []

		for i, variant in enumerate(variants):
			if variant > len(TIME_STEP_FACTOR):
				continue
			x_vals.append(variant)

			time_steps = []
			doubling_times = []
			execution_times = []
			rnap_elong_rates = []
			rnap_fraction_active = []
			rib_elong_rates = []
			rib_fraction_active = []
			protein_fold_changes = []
			dna_fold_changes = []
			sm_fold_changes = []
			trna_fold_changes = []
			rrna_fold_changes = []
			mrna_fold_changes = []
			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				main_reader = TableReader(os.path.join(simOutDir, 'Main'))
				eval_time_reader = TableReader(os.path.join(simOutDir, 'EvaluationTime'))
				ribosome_reader = TableReader(os.path.join(simOutDir, 'RibosomeData'))
				rnap_reader = TableReader(os.path.join(simOutDir, 'RnapData'))
				unique_reader = TableReader(os.path.join(simOutDir, 'UniqueMoleculeCounts'))
				mass_reader = TableReader(os.path.join(simOutDir, 'Mass'))

				# Load data
				time = main_reader.readColumn('time')
				time_step = main_reader.readColumn('timeStepSec')
				clock_time = eval_time_reader.readColumn('clock_time')
				rib_elong_rate = ribosome_reader.readColumn('effectiveElongationRate')
				rnap_elongations = rnap_reader.readColumn('actualElongations')
				unique_ids = unique_reader.readAttribute('uniqueMoleculeIds')
				unique_counts = unique_reader.readColumn('uniqueMoleculeCounts')
				active_rnap_counts = unique_counts[:, unique_ids.index('active_RNAP')]
				active_ribosome_counts = unique_counts[:, unique_ids.index('active_ribosome')]
				protein_mass = mass_reader.readColumn('proteinMass')
				dna_mass = mass_reader.readColumn('dnaMass')
				sm_mass = mass_reader.readColumn('smallMoleculeMass')
				trna_mass = mass_reader.readColumn('tRnaMass')
				rrna_mass = mass_reader.readColumn('rRnaMass')
				mrna_mass = mass_reader.readColumn('mRnaMass')
				(inactive_rnap_counts, inactive_ribosome_counts) = read_bulk_molecule_counts(
					simOutDir, (inactive_rnap_id, ribosome_subunit_ids))

				# Calculate derived values
				rnap_elong_rate = rnap_elongations / active_rnap_counts / time_step
				active_rnap = active_rnap_counts / (active_rnap_counts + inactive_rnap_counts)
				active_ribosomes = active_ribosome_counts / (
					active_ribosome_counts + np.min(inactive_ribosome_counts, axis=1))

				# Save cell data for this variant
				time_steps.append(time_step.mean())
				doubling_times.append((time[-1] - time[0]) / 60)
				execution_times.append((clock_time[-1] - clock_time[0]) / 60)
				rnap_elong_rates.append(rnap_elong_rate.mean())
				rnap_fraction_active.append(active_rnap.mean())
				rib_elong_rates.append(rib_elong_rate.mean())
				rib_fraction_active.append(active_ribosomes.mean())
				protein_fold_changes.append((protein_mass[-1] - protein_mass[0]) / protein_mass[0])
				dna_fold_changes.append((dna_mass[-1] - dna_mass[0]) / dna_mass[0])
				sm_fold_changes.append((sm_mass[-1] - sm_mass[0]) / sm_mass[0])
				trna_fold_changes.append((trna_mass[-1] - trna_mass[0]) / trna_mass[0])
				rrna_fold_changes.append((rrna_mass[-1] - rrna_mass[0]) / rrna_mass[0])
				mrna_fold_changes.append((mrna_mass[-1] - mrna_mass[0]) / mrna_mass[0])

			# Aggregate all data for this variant
			all_time_steps[i] = np.mean(time_steps)
			all_doubling_times[i] = np.mean(doubling_times)
			all_execution_times[i] = np.mean(execution_times)
			all_rnap_elong_rates[i] = np.mean(rnap_elong_rates)
			all_rnap_fraction_active[i] = np.mean(rnap_fraction_active)
			all_rib_elong_rates[i] = np.mean(rib_elong_rates)
			all_rib_fraction_active[i] = np.mean(rib_fraction_active)
			all_protein_fold_changes[i] = np.mean(protein_fold_changes)
			all_dna_fold_changes[i] = np.mean(dna_fold_changes)
			all_sm_fold_changes[i] = np.mean(sm_fold_changes)
			all_trna_fold_changes[i] = np.mean(trna_fold_changes)
			all_rrna_fold_changes[i] = np.mean(rrna_fold_changes)
			all_mrna_fold_changes[i] = np.mean(mrna_fold_changes)

		# Create bar plots
		plt.figure(figsize=(8.5, 11))
		gs = gridspec.GridSpec(5, 3)

		plot_bar(gs[0, 0], x_vals, all_time_steps, 'Time Step\n(sec)')
		plot_bar(gs[0, 1], x_vals, all_doubling_times, 'Doubling Time\n(min)')
		plot_bar(gs[0, 2], x_vals, all_execution_times, 'Execution Time\n(min)')
		plot_bar(gs[1, 0], x_vals, all_protein_fold_changes, 'Protein mass\nfold change')
		plot_bar(gs[1, 1], x_vals, all_dna_fold_changes, 'DNA mass\nfold change')
		plot_bar(gs[1, 2], x_vals, all_sm_fold_changes, 'Small molecule mass\nfold change')
		plot_bar(gs[2, 0], x_vals, all_trna_fold_changes, 'tRNA mass\nfold change')
		plot_bar(gs[2, 1], x_vals, all_rrna_fold_changes, 'rRNA mass\nfold change')
		plot_bar(gs[2, 2], x_vals, all_mrna_fold_changes, 'mRNA mass\nfold change', show_x=True)
		plot_bar(gs[3, 0], x_vals, all_rnap_elong_rates, 'RNAP elongation rate\n(nt/s)')
		plot_bar(gs[3, 1], x_vals, all_rnap_fraction_active, 'RNAP fraction active')
		plot_bar(gs[4, 0], x_vals, all_rib_elong_rates, 'Ribosome elongation rate\n(AA/s)', show_x=True)
		plot_bar(gs[4, 1], x_vals, all_rib_fraction_active, 'Ribosome fraction active', show_x=True)

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
