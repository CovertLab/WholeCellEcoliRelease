"""
Compare various cell macromolecular properties across different growth rates
with validation data from Dennis & Bremmer (2021)
https://journals.asm.org/doi/epub/10.1128/ecosal.5.2.3
"""

import os

import numpy as np
from matplotlib import pyplot as plt
import _pickle as cPickle

from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.io.tablereader import TableReader
from wholecell.utils import units

from models.ecoli.analysis import variantAnalysisPlot
from wholecell.analysis.analysis_tools import (exportFigure, read_stacked_columns)

FIG_SIZE = (10, 5)

LABEL_FONT_SIZE=9
LEGEND_FONT_SIZE = 6

SIM_LINE_STYLE = dict(
    color = "black",
    fmt = ' ',
    marker = 'o',
    markersize = 1,
    linewidth = 0.5,
    linestyle = 'dashed'
)

VAL_LINE_STYLE = dict(
    color = "blue",
    marker = 'o',
    markersize = 1,
    linewidth = 0.5,
    linestyle = 'dashed'
)

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
    def plot_lines(self, ax, sim_x, sim_y, sim_std, val_x, val_y, title, y_lim=None):
        # Sort plotted variables in order of increasing doubling time
        sim_y = sim_y[np.argsort(sim_x)][::-1]
        sim_std = sim_std[np.argsort(sim_x)][::-1]
        sim_x = sim_x[np.argsort(sim_x)][::-1]
        val_y = val_y[np.argsort(val_x)][::-1]
        val_x = val_x[np.argsort(val_x)][::-1]

        # Make plot
        ax.errorbar(sim_x, sim_y, yerr=sim_std, label = "Simulation", **SIM_LINE_STYLE)
        ax.errorbar(val_x, val_y, yerr=0, label = "Bremer & Dennis 2021", **VAL_LINE_STYLE)
        if y_lim is not None:
            ax.set_ylim(y_lim)
        ax.set_xlabel("Doubling time (min)", fontsize=LABEL_FONT_SIZE)
        ax.set_title(title, fontsize=LABEL_FONT_SIZE)
        ax.legend(fontsize=LEGEND_FONT_SIZE)                    

    def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        ap = AnalysisPaths(inputDir, variant_plot=True)

        # Plotted units of data
        mass_units = units.fg
        time_units = units.min

        # Extract validation data and strip the units while numerically
        # converting to the plotted units (validationDataFile stores a series of
        # ndarrays, each containing the data for a particular value in
        # increasing order of growth rate, and each ndarray tagged with its
        # corresponding unit)
        validation_data = cPickle.load(open(validationDataFile, "rb"))
        db_table = validation_data.macromolecular_growth_rate_modulation
        val_doubling_time = db_table.doubling_time.asNumber(time_units)
        val_PRD_per_mass = db_table.PRD_per_mass.asNumber()
        val_dry_mass = db_table.mass_per_cell.asNumber(mass_units)
        val_protein_mass_per_cell = db_table.protein_per_cell_ug.asNumber(mass_units)
        val_RNA_mass_per_cell = db_table.RNA_per_cell_ug.asNumber(mass_units)
        val_DNA_mass_per_cell = db_table.DNA_per_cell_ug.asNumber(mass_units)
        val_PRD_per_cell = db_table.PRD_per_cell.asNumber(mass_units)
        val_total_RNA_stable_fraction = db_table.total_RNA_stable_fraction.asNumber()
        val_stable_RNA_tRNA_fraction = db_table.stable_RNA_tRNA_fraction.asNumber()

        # Initialize ndarrays for simulation data means and standard deviations
        # across the cells of each variant
        sim_doubling_time = np.zeros(ap.n_variant)
        sim_PRD_per_mass = np.zeros(ap.n_variant)
        sim_dry_mass = np.zeros(ap.n_variant)
        sim_protein_mass = np.zeros(ap.n_variant)
        sim_RNA_mass = np.zeros(ap.n_variant)
        sim_DNA_mass = np.zeros(ap.n_variant)
        sim_PRD_per_cell = np.zeros(ap.n_variant)
        sim_total_RNA_stable_fraction = np.zeros(ap.n_variant)
        sim_stable_RNA_tRNA_fraction = np.zeros(ap.n_variant)

        sim_PRD_per_mass_std = np.zeros(ap.n_variant)
        sim_dry_mass_std = np.zeros(ap.n_variant)
        sim_protein_mass_std = np.zeros(ap.n_variant)
        sim_RNA_mass_std = np.zeros(ap.n_variant)
        sim_DNA_mass_std = np.zeros(ap.n_variant)
        sim_PRD_per_cell_std = np.zeros(ap.n_variant)
        sim_total_RNA_stable_fraction_std = np.zeros(ap.n_variant)
        sim_stable_RNA_tRNA_fraction_std = np.zeros(ap.n_variant)

        variants = ap.get_variants()

        # Read table attributes such as units from a test cell for later use
        test_cell = ap.get_cells(variant=[variants[0]])[0]
        simOutDir = os.path.join(test_cell, "simOut")
        try:
            mass = TableReader(os.path.join(simOutDir, "Mass"))
        except Exception as e:
            print('Error with data for %s: %s' % (test_cell, e))
        protein_mass_units = getattr(units, mass.readAttribute("protein_units"))
        nucleic_acid_mass_units = getattr(units, mass.readAttribute("rna_units"))
        dry_mass_units = getattr(units, mass.readAttribute("cellDry_units"))

        for varIdx in range(ap.n_variant):
            variant = variants[varIdx]
            cells = ap.get_cells(variant=[variant])

            # Read data from the cell's listeners
            doubling_time = read_stacked_columns(cells, 'Main', 'time', fun=lambda x: x[-1] - x[0]) * (units.s).asNumber(time_units)
            protein_mass = read_stacked_columns(cells, 'Mass', 'proteinMass', fun=lambda x: x.mean()) * protein_mass_units.asNumber(mass_units)
            mRna_mass = read_stacked_columns(cells, 'Mass', 'mRnaMass', fun=lambda x: x.mean()) * nucleic_acid_mass_units.asNumber(mass_units)
            tRna_mass = read_stacked_columns(cells, 'Mass', 'tRnaMass', fun=lambda x: x.mean()) * nucleic_acid_mass_units.asNumber(mass_units)
            rRna_mass = read_stacked_columns(cells, 'Mass', 'rRnaMass', fun=lambda x: x.mean()) * nucleic_acid_mass_units.asNumber(mass_units)
            rna_mass = read_stacked_columns(cells, 'Mass', 'rnaMass', fun=lambda x: x.mean()) * nucleic_acid_mass_units.asNumber(mass_units)
            dna_mass = read_stacked_columns(cells, 'Mass', 'dnaMass', fun=lambda x: x.mean()) * nucleic_acid_mass_units.asNumber(mass_units)
            dry_mass = read_stacked_columns(cells, 'Mass', 'dryMass', fun=lambda x: x.mean()) * dry_mass_units.asNumber(mass_units)
            PRD_per_mass = (protein_mass + rna_mass + dna_mass) / dry_mass
            PRD_per_cell = protein_mass + rna_mass + dna_mass
            # In Dennis & Bremmer (2021), the stable fraction of total RNA
            # is calculated by estimating the fraction of total RNA that is
            # mRNA, then subtracting from 1. In the model, unstable RNAs
            # also include miscellaneous RNAs not accounted for in the mRNA
            # count. We thus use the more "accurate" measure of tRNA + rRNA,
            # though the difference between the two methods is quite small,
            # at most 0.01.
            total_rna_stable_fraction = (tRna_mass + rRna_mass) / rna_mass
            stable_rna_tRna_fraction = tRna_mass / (tRna_mass + rRna_mass)

            # Calculate mean and standard deviation for each value across
            # the cells of the current variant
            sim_doubling_time[varIdx] = np.mean(doubling_time)
            sim_PRD_per_mass[varIdx] = np.mean(PRD_per_mass)
            sim_dry_mass[varIdx] = np.mean(dry_mass)
            sim_protein_mass[varIdx] = np.mean(protein_mass)
            sim_RNA_mass[varIdx] = np.mean(rna_mass)
            sim_DNA_mass[varIdx] = np.mean(dna_mass)
            sim_PRD_per_cell[varIdx] = np.mean(PRD_per_cell)
            sim_total_RNA_stable_fraction[varIdx] = np.mean(total_rna_stable_fraction)
            sim_stable_RNA_tRNA_fraction[varIdx] = np.mean(stable_rna_tRna_fraction)

            sim_PRD_per_mass_std[varIdx] = PRD_per_mass.std()
            sim_dry_mass_std[varIdx] = dry_mass.std()
            sim_protein_mass_std[varIdx] = protein_mass.std()
            sim_RNA_mass_std[varIdx] = rna_mass.std()
            sim_DNA_mass_std[varIdx] = dna_mass.std()
            sim_PRD_per_cell_std[varIdx] = PRD_per_cell.std()
            sim_total_RNA_stable_fraction_std[varIdx] = total_rna_stable_fraction.std()
            sim_stable_RNA_tRNA_fraction_std[varIdx] = stable_rna_tRna_fraction.std()

        # Make plots
        fig, axes  = plt.subplots(ncols=4, nrows=2, figsize=FIG_SIZE)

        self.plot_lines(axes[0, 0], sim_doubling_time, sim_PRD_per_mass, sim_PRD_per_mass_std,
                        val_doubling_time, val_PRD_per_mass, "Protein + RNA + DNA mass per dry mass", y_lim=[0,1])
        self.plot_lines(axes[0, 1], sim_doubling_time, sim_dry_mass, sim_dry_mass_std,
                        val_doubling_time, val_dry_mass, "Dry mass per cell (fg)")
        self.plot_lines(axes[0, 2], sim_doubling_time, sim_total_RNA_stable_fraction, sim_total_RNA_stable_fraction_std,
                        val_doubling_time, val_total_RNA_stable_fraction, "Stable Fraction of Total RNA", y_lim=[0.8, 1])
        self.plot_lines(axes[0, 3], sim_doubling_time, sim_stable_RNA_tRNA_fraction, sim_stable_RNA_tRNA_fraction_std,
                        val_doubling_time, val_stable_RNA_tRNA_fraction, "tRNA Fraction of Stable RNA", y_lim=[0, 0.2])
        self.plot_lines(axes[1, 0], sim_doubling_time, sim_protein_mass, sim_protein_mass_std,
                        val_doubling_time, val_protein_mass_per_cell, "Protein mass per cell (fg)")
        self.plot_lines(axes[1, 1], sim_doubling_time, sim_RNA_mass, sim_RNA_mass_std,
                        val_doubling_time, val_RNA_mass_per_cell, "RNA mass per cell (fg)")
        self.plot_lines(axes[1, 2], sim_doubling_time, sim_DNA_mass, sim_DNA_mass_std,
                        val_doubling_time, val_DNA_mass_per_cell, "DNA mass per cell (fg)")
        self.plot_lines(axes[1, 3], sim_doubling_time, sim_PRD_per_cell, sim_PRD_per_cell_std,
                        val_doubling_time, val_PRD_per_cell, "Protein + RNA + DNA mass per cell (fg)")

        plt.tight_layout()

        exportFigure(plt, plotOutDir, plotOutFileName, metadata)

if __name__ == "__main__":
    Plot().cli()
