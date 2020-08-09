from __future__ import absolute_import


import os
import pickle

import numpy as np
from matplotlib import pyplot as plt

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader

from wholecell.utils.sparkline import whitePadSparklineAxis
from wholecell.analysis.analysis_tools import exportFigure
from models.ecoli.analysis import variantAnalysisPlot

FONT_SIZE=9

SIM_PLOT_STYLE = dict(
    label="Simulation",
    fmt='',
    marker='o',
    markersize=5,
    linewidth=1,
    capsize=2)

EXP_PLOT_STYLE = dict(
    label="Experimental",
    marker='o',
    markersize=5,
    linewidth=1,
    )

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
    def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        if not os.path.isdir(inputDir):
            raise Exception, "variantDir does not currently exist as a directory"

        if not os.path.exists(plotOutDir):
            os.mkdir(plotOutDir)

        ap = AnalysisPaths(inputDir, variant_plot = True)

        fig = plt.figure()
        fig.set_figwidth(5)
        fig.set_figheight(5)

        bremer_tau = [40, 100, 24]

        bremer_origins_per_cell_at_initiation = [2, 1, 4]
        bremer_rrn_init_rate = [20*23, 4*12.4, 58*35.9]

        bremer_rna_mass_per_cell = [77, 20,  211]
        bremer_elng_rate = [18, 12,  21]

        sim_doubling_time = np.zeros(ap.n_variant)
        sim_doubling_time_std = np.zeros(ap.n_variant)

        sim_origins_per_cell_at_initiation = np.zeros(ap.n_variant)
        sim_rna_mass_per_cell = np.zeros(ap.n_variant)
        sim_elng_rate = np.zeros(ap.n_variant)
        sim_rrn_init_rate = np.zeros(ap.n_variant)

        sim_origins_per_cell_at_initiation_std = np.zeros(ap.n_variant)
        sim_elng_rate_std = np.zeros(ap.n_variant)
        sim_rna_mass_per_cell_std = np.zeros(ap.n_variant)
        sim_rrn_init_rate_std = np.zeros(ap.n_variant)

        variants = ap.get_variants()

        for varIdx in range(ap.n_variant):
            variant = variants[varIdx]

            print("variant {}".format(variant))

            all_cells = ap.get_cells(variant=[variant])

            print("Total cells: {}".format(len(all_cells)))

            try:
                sim_data = pickle.load(open(ap.get_variant_kb(variant)))
            except Exception as e:
                print "Couldn't load sim_data object. Exiting.", e
                return

            num_origin_at_init = np.zeros(len(all_cells))
            doubling_time = np.zeros(len(all_cells))
            meanRnaMass = np.zeros(len(all_cells))
            meanElngRate = np.zeros(len(all_cells))
            meanRrnInitRate = np.zeros(len(all_cells))

            for idx, simDir in enumerate(all_cells):
                print "cell {} of {}".format(idx, len(all_cells))

                simOutDir = os.path.join(simDir, "simOut")

                try:
                    time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")
                    doubling_time[idx] = time[-1] - time[0]
                except Exception as e:
                    print 'Error with data for %s: %s' % (simDir, e)
                    continue

                timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")

                meanRnaMass[idx] = TableReader(os.path.join(simOutDir, "Mass")).readColumn("rnaMass").mean()
                meanElngRate[idx] = TableReader(os.path.join(simOutDir, "RibosomeData")).readColumn("effectiveElongationRate").mean()

                numOrigin = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("numberOfOric")

                massPerOric = TableReader(os.path.join(simOutDir, "ReplicationData")).readColumn("criticalMassPerOriC")
                idxInit = np.where(massPerOric >= 1)[0]
                numOriginAtInit = numOrigin[idxInit - 1]
                if numOriginAtInit.size:
                    num_origin_at_init[idx] = numOriginAtInit.mean()
                else:
                    num_origin_at_init[idx] = np.nan

                transcriptDataFile = TableReader(os.path.join(simOutDir, "TranscriptElongationListener"))
                rnaSynth = transcriptDataFile.readColumn("countRnaSynthesized")
                isRRna = sim_data.process.transcription.rnaData["isRRna"]
                meanRrnInitRate[idx] = (rnaSynth[:, isRRna].sum(axis=1) / timeStepSec).mean() * 60. / 3

            sim_rna_mass_per_cell[varIdx] = meanRnaMass.mean()
            sim_elng_rate[varIdx] = meanElngRate.mean()
            sim_origins_per_cell_at_initiation[varIdx] = np.nanmean(num_origin_at_init)
            sim_doubling_time[varIdx] = np.nanmean(doubling_time) / 60.
            sim_rrn_init_rate[varIdx] = np.nanmean(meanRrnInitRate)

            sim_rna_mass_per_cell_std[varIdx] = meanRnaMass.std()
            sim_elng_rate_std[varIdx] = meanElngRate.std()
            sim_origins_per_cell_at_initiation_std[varIdx] = np.nanstd(num_origin_at_init)
            sim_doubling_time_std[varIdx] = np.nanstd(doubling_time) / 60.
            sim_rrn_init_rate_std[varIdx] = np.nanstd(meanRrnInitRate)

        bremer_tau = np.array(bremer_tau)

        ax0 = plt.subplot2grid((2,2), (0,0))
        ax1 = plt.subplot2grid((2,2), (1,0), sharex=ax0)
        ax2 = plt.subplot2grid((2,2), (0,1), sharex=ax0)
        ax3 = plt.subplot2grid((2,2), (1,1), sharex=ax0)

        lines = {'linestyle': 'dashed'}
        plt.rc('lines', **lines)
        plt.style.use('seaborn-deep')
        color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

        ax0.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]], sim_rna_mass_per_cell[np.argsort(sim_doubling_time)[::-1]], yerr=sim_rna_mass_per_cell_std[np.argsort(sim_doubling_time)[::-1]], color=color_cycle[0], **SIM_PLOT_STYLE)
        ax0.errorbar(bremer_tau[np.argsort(bremer_tau)[::-1]], np.array(bremer_rna_mass_per_cell)[np.argsort(bremer_tau)[::-1]], color=color_cycle[2], **EXP_PLOT_STYLE)
        ax0.set_title("RNA mass per cell (fg)", fontsize=FONT_SIZE)
        ax0.set_xlim([0, 135])
        ax0.set_ylim([0, 250])
        ax0.legend(loc=1, fontsize='xx-small', markerscale=0.5, frameon=False)

        ax1.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]], sim_elng_rate[np.argsort(sim_doubling_time)[::-1]], yerr=sim_elng_rate_std[np.argsort(sim_doubling_time)[::-1]], color=color_cycle[0], **SIM_PLOT_STYLE)
        ax1.errorbar(bremer_tau[np.argsort(bremer_tau)[::-1]], np.array(bremer_elng_rate)[np.argsort(bremer_tau)[::-1]], color=color_cycle[2], **EXP_PLOT_STYLE)
        ax1.set_title("Ribosome elongation\nrate (aa/s/ribosome)", fontsize=FONT_SIZE)
        ax1.set_xlabel("Doubling time (min)", fontsize=FONT_SIZE)
        ax1.set_ylim([0, 24])

        ax2.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]], sim_origins_per_cell_at_initiation[np.argsort(sim_doubling_time)[::-1]], yerr=sim_origins_per_cell_at_initiation_std[np.argsort(sim_doubling_time)[::-1]], color=color_cycle[0], **SIM_PLOT_STYLE)
        ax2.errorbar(bremer_tau[np.argsort(bremer_tau)[::-1]], np.array(bremer_origins_per_cell_at_initiation)[np.argsort(bremer_tau)[::-1]], color=color_cycle[2], **EXP_PLOT_STYLE)
        ax2.set_title("Average origins at chrom. init.", fontsize=FONT_SIZE)
        ax2.set_ylim([0.5, 4.5])

        ax3.errorbar(sim_doubling_time[np.argsort(sim_doubling_time)[::-1]], sim_rrn_init_rate[np.argsort(sim_doubling_time)[::-1]], yerr=sim_rrn_init_rate_std[np.argsort(sim_doubling_time)[::-1]], color=color_cycle[0], **SIM_PLOT_STYLE)
        ax3.errorbar(bremer_tau[np.argsort(bremer_tau)[::-1]], np.array(bremer_rrn_init_rate)[np.argsort(bremer_tau)[::-1]], color=color_cycle[2], **EXP_PLOT_STYLE)
        ax3.set_title("Rate of rrn initiation (1/min)", fontsize=FONT_SIZE)
        ax3.set_ylim([0, 2500])

        # ax3.legend(loc=1, frameon=True, fontsize=7)
        ax3.set_xlabel("Doubling time (min)", fontsize=FONT_SIZE)

        axes_list = [ax0, ax1, ax2, ax3]

        for a in axes_list:
            for tick in a.yaxis.get_major_ticks():
                tick.label.set_fontsize(FONT_SIZE)
            for tick in a.xaxis.get_major_ticks():
                tick.label.set_fontsize(FONT_SIZE)

        whitePadSparklineAxis(ax0, False)
        whitePadSparklineAxis(ax1)
        whitePadSparklineAxis(ax2, False)
        whitePadSparklineAxis(ax3)

        plt.subplots_adjust(bottom = 0.2, wspace=0.3)

        exportFigure(plt, plotOutDir, plotOutFileName, metadata)


if __name__ == "__main__":
    Plot().cli()
