from __future__ import absolute_import

import os
import cPickle
import numpy as np
from matplotlib import pyplot as plt
from multiprocessing import Pool

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.utils import parallelization
from models.ecoli.analysis import variantAnalysisPlot

FONT_SIZE =9
X_LIM = [20, 120]
HIGHLIGHT_COLOR = 'tab:orange'

SIM_PLOT_STYLE = dict(
	label='Simulation',
	fmt='',
	marker='o',
	markersize=10,
	linewidth=1,
	capsize=4,
	elinewidth=2)

EXP_PLOT_STYLE = dict(
	label='Experiment',
	marker='o',
	markersize=10,
	linewidth=1,
	)

# Measurements from Bremer
bremer_tau = np.array([40, 100, 24])
bremer_origins_per_cell_at_initiation = np.array([2, 1, 4])
bremer_rrn_init_rate = np.array([20*23, 4*12.4, 58*35.9])
bremer_rna_mass_per_cell = np.array([77, 20,  211])
bremer_elng_rate = np.array([18, 12,  21])


def mp_worker(sim_dir):
	sim_out_dir = os.path.join(sim_dir, 'simOut')
	doubling_time = np.nan
	mean_rna_mass = np.nan
	mean_elng_rate = np.nan
	n_origin_init = np.nan
	mean_rrn_init_rate = np.nan

	try:
		# Doubling time
		time = TableReader(os.path.join(sim_out_dir, "Main")).readColumn("time")
		doubling_time = time[-1] - time[0]

		# RNA mass
		mean_rna_mass = TableReader(os.path.join(sim_out_dir, "Mass")).readColumn("rnaMass").mean()

		# Ribosome elongation rate
		mean_elng_rate = TableReader(
			os.path.join(sim_out_dir, "RibosomeData")).readColumn("effectiveElongationRate").mean()

		# Number of origins at chromosome initiation
		n_origin = TableReader(os.path.join(sim_out_dir, "ReplicationData")).readColumn("numberOfOric")
		mass_per_oric = TableReader(
			os.path.join(sim_out_dir, "ReplicationData")).readColumn("criticalMassPerOriC")
		index_init = np.where(mass_per_oric >= 1)[0]
		n_origin_init_raw = n_origin[index_init - 1]
		if n_origin_init_raw.size:
			n_origin_init = n_origin_init_raw.mean()

		# rRNA initiation rate
		rnaSynth = TableReader(
			os.path.join(sim_out_dir, "TranscriptElongationListener")).readColumn("countRnaSynthesized")
		time_step_sec = TableReader(os.path.join(sim_out_dir, "Main")).readColumn("timeStepSec")
		mean_rrn_init_rate = (rnaSynth[:, is_rRNA].sum(axis=1) / time_step_sec).mean() * 60. / 3

	except Exception as e:
		print('Excluded from analysis due to broken files: {}'.format(sim_out_dir))

	return [doubling_time, mean_rna_mass, mean_elng_rate, n_origin_init, mean_rrn_init_rate]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, "variantDir does not currently exist as a directory"

		if not os.path.exists(plotOutDir):
			os.mkdir(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot = True)
		variants = ap.get_variants()

		index_doubling_time = 0
		sim_doubling_time = []

		index_rna_mass = 1
		sim_rna_mass_per_cell = []
		sim_rna_mass_per_cell_std = []

		index_elng_rate = 2
		sim_elng_rate = []
		sim_elng_rate_std = []

		index_n_origin_init = 3
		sim_origins_per_cell_at_initiation = []
		sim_origins_per_cell_at_initiation_std = []

		index_rrn_init_rate = 4
		sim_rrn_init_rate = []
		sim_rrn_init_rate_std = []

		for varIdx in range(ap.n_variant):
			variant = variants[varIdx]
			print("variant {}".format(variant))

			sim_dirs = ap.get_cells(variant=[variant])
			n_sims = len(sim_dirs)
			print("Total cells: {}".format(n_sims))

			try:
				sim_data = cPickle.load(open(ap.get_variant_kb(variant)))

				global is_rRNA
				is_rRNA = sim_data.process.transcription.rnaData["isRRna"]

			except Exception as e:
				print "Couldn't load sim_data object. Exiting.", e
				return

			p = Pool(parallelization.cpus())
			output = np.array(p.map(mp_worker, sim_dirs))
			p.close()
			p.join()

			# Filter output from broken files using np.nanmean and np.nanstd
			sim_doubling_time.append(np.nanmean(output[:, index_doubling_time]) / 60.)

			sim_rna_mass_per_cell.append(np.nanmean(output[:, index_rna_mass]))
			sim_rna_mass_per_cell_std.append(np.nanstd(output[:, index_rna_mass]))

			sim_elng_rate.append(np.nanmean(output[:, index_elng_rate]))
			sim_elng_rate_std.append(np.nanstd(output[:, index_elng_rate]))

			sim_origins_per_cell_at_initiation.append(np.nanmean(output[:, index_n_origin_init]))
			sim_origins_per_cell_at_initiation_std.append(np.nanstd(output[:, index_n_origin_init]))

			sim_rrn_init_rate.append(np.nanmean(output[:, index_rrn_init_rate]))
			sim_rrn_init_rate_std.append(np.nanstd(output[:, index_rrn_init_rate]))

		sim_doubling_time = np.array(sim_doubling_time)

		# Plot
		fig, axes_list = plt.subplots(1, 4, figsize=(15, 5))
		ax0, ax1, ax2, ax3 = axes_list
		sort_sim = np.argsort(sim_doubling_time)[::-1]
		sort_bremer = np.argsort(bremer_tau)[::-1]

		# RNA mass per cell
		ax0.errorbar(
			sim_doubling_time[sort_sim],
			np.array(sim_rna_mass_per_cell)[sort_sim],
			yerr=np.array(sim_rna_mass_per_cell_std)[sort_sim],
			color='tab:blue', **SIM_PLOT_STYLE)
		ax0.errorbar(
			bremer_tau[sort_bremer],
			bremer_rna_mass_per_cell[sort_bremer],
			color=HIGHLIGHT_COLOR, **EXP_PLOT_STYLE)
		ax0.set_title('RNA mass per cell (fg)', fontsize=FONT_SIZE)
		ax0.set_xlabel('Doubling time (min)', fontsize=FONT_SIZE)
		ax0.set_xlim([0, 135])
		ax0.set_ylim([0, 250])
		ax0.legend(loc=1, fontsize='xx-small', markerscale=0.5, frameon=False)

		# Ribosome elongation rate
		ax1.errorbar(
			sim_doubling_time[sort_sim],
			np.array(sim_elng_rate)[sort_sim],
			yerr=np.array(sim_elng_rate_std)[sort_sim],
			color='tab:blue', **SIM_PLOT_STYLE)
		ax1.errorbar(
			bremer_tau[sort_bremer],
			bremer_elng_rate[sort_bremer],
			color=HIGHLIGHT_COLOR, **EXP_PLOT_STYLE)
		ax1.set_title('Ribosome elongation\nrate (aa/s/ribosome)', fontsize=FONT_SIZE)
		ax1.set_xlabel('Doubling time (min)', fontsize=FONT_SIZE)
		ax1.set_ylim([5, 24])

		# Number of origins at chromosome initiation
		ax2.errorbar(
			sim_doubling_time[sort_sim],
			np.array(sim_origins_per_cell_at_initiation)[sort_sim],
			yerr=np.array(sim_origins_per_cell_at_initiation_std)[sort_sim],
			color='tab:blue', **SIM_PLOT_STYLE)
		ax2.errorbar(
			bremer_tau[sort_bremer],
			bremer_origins_per_cell_at_initiation[sort_bremer],
			color=HIGHLIGHT_COLOR, **EXP_PLOT_STYLE)
		ax2.set_title('Average origins at chrom. init.', fontsize=FONT_SIZE)
		ax2.set_xlabel('Doubling time (min)', fontsize=FONT_SIZE)
		ax2.set_ylim([0.5, 4.5])

		# rRNA initiation rate
		ax3.errorbar(
			sim_doubling_time[sort_sim],
			np.array(sim_rrn_init_rate)[sort_sim],
			yerr=np.array(sim_rrn_init_rate_std)[sort_sim],
			color='tab:blue', **SIM_PLOT_STYLE)
		ax3.errorbar(
			bremer_tau[sort_bremer],
			bremer_rrn_init_rate[sort_bremer],
			color=HIGHLIGHT_COLOR, **EXP_PLOT_STYLE)
		ax3.set_title('Rate of rrn initiation (1/min)', fontsize=FONT_SIZE)
		ax3.set_ylim([0, 2500])
		ax3.set_xlabel('Doubling time (min)', fontsize=FONT_SIZE)

		for ax in axes_list:
			ax.set_xlim(X_LIM)
			ax.set_xticks(X_LIM)
			ax.set_ylim(ax.get_ylim())
			ax.set_yticks(ax.get_ylim())

			for tick in ax.yaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)
			for tick in ax.xaxis.get_major_ticks():
				tick.label.set_fontsize(FONT_SIZE)

		plt.subplots_adjust(bottom=0.25, top=0.75, left=0.05, right=0.95, wspace=0.4)
		exportFigure(plt, plotOutDir, '{}__test'.format(plotOutFileName), metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
