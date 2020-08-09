"""
Analysis of fluxes for the factorial design experiments generated by the
variant kinetic_constraints_factorial_experiments.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import pickle
import os

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import pearsonr

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS, MASS_UNITS
from models.ecoli.sim.variants.kinetic_constraints_factorial_experiments import get_disabled_constraints
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, filepath, units


# IDs of interest
SUCC_ID = 'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.'
NADH_ID = 'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)'
GLC_ID = 'GLC[p]'
REACTIONS = [SUCC_ID, NADH_ID]
EXCHANGES = [GLC_ID]
HIGHLIGHTED_CONSTRAINTS = [
    SUCC_ID,
    NADH_ID,
    'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.',
    'GLUTATHIONE-REDUCT-NADPH-RXN',
    ]

# Flux units
MODEL_FLUX_UNITS = COUNTS_UNITS / MASS_UNITS / TIME_UNITS
OUTPUT_FLUX_UNITS = units.mmol / units.g / units.h
FLUX_CONVERSION = MODEL_FLUX_UNITS.asNumber(OUTPUT_FLUX_UNITS)

# Plot region cutoffs
GLC_MAX = 11.5
SUCC_DISTANCE = np.log2(2)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
    def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
        if not os.path.isdir(inputDir):
            raise Exception, 'inputDir does not currently exist as a directory'

        filepath.makedirs(plotOutDir)

        ap = AnalysisPaths(inputDir, variant_plot=True)
        variants = ap.get_variants()
        n_variants = len(variants)

        # Load sim_data
        with open(os.path.join(inputDir, 'kb', constants.SERIALIZED_FIT1_FILENAME), 'rb') as f:
            sim_data = pickle.load(f)
        cell_density = sim_data.constants.cellDensity.asNumber(MASS_UNITS / VOLUME_UNITS)

        # Load validation_data
        with open(validationDataFile, "rb") as f:
            validation_data = pickle.load(f)
        toyaReactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
        toyaFluxes = validation_data.reactionFlux.toya2010fluxes["reactionFlux"]
        toyaStdev = validation_data.reactionFlux.toya2010fluxes["reactionFluxStdev"]
        toyaFluxesDict = dict(zip(toyaReactions, toyaFluxes))
        toyaStdevDict = dict(zip(toyaReactions, toyaStdev))

        glc_uptakes = np.zeros(n_variants)
        log_ratio_succ = np.zeros(n_variants)
        size_pearson = np.zeros(n_variants)
        selected_indicies = np.zeros(n_variants, bool)
        for v, variant in enumerate(variants):
            # initialize kinetic flux comparison
            exchange_fluxes = {entry: [] for entry in EXCHANGES}
            reaction_fluxes = {entry: [] for entry in REACTIONS}

            modelFluxes = {}
            toyaOrder = []
            for rxn in toyaReactions:
                modelFluxes[rxn] = []
                toyaOrder.append(rxn)

            for sim_dir in ap.get_cells(variant=[variant]):
                simOutDir = os.path.join(sim_dir, "simOut")

                try:
                    # Listeners used
                    massListener = TableReader(os.path.join(simOutDir, "Mass"))
                    fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
                    enzymeKineticsReader = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))

                    ## Read from mass listener
                    cellMass = massListener.readColumn("cellMass")
                    # skip if no data
                    if cellMass.shape is ():
                        continue
                    dryMass = massListener.readColumn("dryMass")
                except Exception as e:
                    print(e)
                    continue

                coefficient = (dryMass / cellMass * cell_density).reshape(-1, 1)

                ## Read from FBA listener
                reactionIDs = {r: i for i, r in enumerate(fbaResults.readAttribute("reactionIDs"))}
                exMolec = {m: i for i, m in enumerate(fbaResults.readAttribute("externalMoleculeIDs"))}
                reactionFluxes = FLUX_CONVERSION * (fbaResults.readColumn("reactionFluxes") / coefficient)[1:, :]
                exFlux = fbaResults.readColumn("externalExchangeFluxes")[1:, :]

                ## Read from EnzymeKinetics listener
                constrainedReactions = {r: i for i, r in enumerate(enzymeKineticsReader.readAttribute("constrainedReactions"))}

                ## Append values for relevant reactions.
                # append to exchanges
                for entry in EXCHANGES:
                    exchange_fluxes[entry].extend(list(exFlux[:, exMolec[entry]]))
                # append to reaction fluxes
                for entry in REACTIONS:
                    reaction_fluxes[entry].extend(list(reactionFluxes[:, reactionIDs[entry]]))

                ## get all Toya reactions, and corresponding simulated fluxes.
                toya_idx = {r: [] for r in toyaReactions}
                for rxn, i in reactionIDs.items():
                    rxn = rxn.split(' (reverse)')
                    if len(rxn) > 1:
                        i = -i
                    rxn = rxn[0].split('__')[0]
                    if rxn in toya_idx:
                        toya_idx[rxn] += [i]
                for toyaReaction, reaction_idx in toya_idx.items():
                    flux_time_course = np.sum([np.sign(i) * reactionFluxes[:, np.abs(i)] for i in reaction_idx], axis=0)
                    modelFluxes[toyaReaction].append(flux_time_course.mean())

            ## Flux comparison with Toya
            toyaVsReactionAve = []
            rxn_order = []
            for rxn, toyaFlux in toyaFluxesDict.iteritems():
                rxn_order.append(rxn)
                if rxn in modelFluxes:
                    toyaVsReactionAve.append((np.mean(modelFluxes[rxn]), toyaFlux.asNumber(OUTPUT_FLUX_UNITS), np.std(modelFluxes[rxn]), toyaStdevDict[rxn].asNumber(OUTPUT_FLUX_UNITS)))

            toyaVsReactionAve = np.array(toyaVsReactionAve)
            rWithAll = pearsonr(toyaVsReactionAve[:,0], toyaVsReactionAve[:,1])
            succ_toya_flux = toyaVsReactionAve[rxn_order.index(SUCC_ID), 1]

            # Save data for plotting
            glc_uptakes[v] = -np.mean(exchange_fluxes[GLC_ID])
            log_ratio_succ[v] = np.log2(np.mean(reaction_fluxes[SUCC_ID]) / succ_toya_flux)
            size_pearson[v] = (rWithAll[0] * 8)**2
            selected_indicies[v] = np.all([c not in constrainedReactions for c in HIGHLIGHTED_CONSTRAINTS])

        # Plot scatterplot
        fig = plt.figure(figsize=(5, 5))
        gs = gridspec.GridSpec(40, 40)

        ## Plot full data
        plt.scatter(glc_uptakes[~selected_indicies], log_ratio_succ[~selected_indicies],
            color='blue', alpha=0.6, s=size_pearson[~selected_indicies])
        plt.scatter(glc_uptakes[selected_indicies], log_ratio_succ[selected_indicies],
            color='red', alpha=0.6, s=size_pearson[selected_indicies])
        x_min, x_max = plt.xlim()
        y_max = max(np.abs(plt.ylim()))
        plt.axvspan(0, GLC_MAX, facecolor='g', alpha=0.1)
        plt.axhspan(-SUCC_DISTANCE, SUCC_DISTANCE, facecolor='g', alpha=0.1)
        plt.axhline(y=0, color='k', linestyle='--')

        ## Format axes
        plt.ylabel('log2(model flux / Toya flux)')
        plt.xlabel('glucose uptake (mmol / g DCW / hr)')
        plt.xlim([np.floor(min(x_min, 10)), np.ceil(x_max)])
        plt.ylim([-y_max, y_max])

        ## Plot highlighted region data
        fig.add_subplot(gs[1:28, -20:-1])
        in_region = (glc_uptakes < GLC_MAX) & (np.abs(log_ratio_succ) < SUCC_DISTANCE)
        selected_in = in_region & selected_indicies
        not_selected_in = in_region & ~selected_indicies
        constraint_labels = np.array([
            [c[:2] for c in constraints] if constraints is not None else []
            for _, constraints in map(get_disabled_constraints, variants)
            ])
        plt.scatter(glc_uptakes[not_selected_in], log_ratio_succ[not_selected_in],
            color='blue', alpha=0.6, s=size_pearson[not_selected_in])
        plt.scatter(glc_uptakes[selected_in], log_ratio_succ[selected_in],
            color='red', alpha=0.6, s=size_pearson[selected_in])
        for x, y, label in zip(glc_uptakes[in_region], log_ratio_succ[in_region], constraint_labels[in_region]):
            plt.text(x, y, ', '.join(label), ha='center', va='top', fontsize=6)
        x_min, _ = plt.xlim()
        x_min = np.floor(min(x_min, 10))
        plt.axvspan(x_min, GLC_MAX, facecolor='g', alpha=0.1)
        plt.axhspan(-SUCC_DISTANCE, SUCC_DISTANCE, facecolor='g', alpha=0.1)

        ## Format axes
        plt.xlim([x_min, GLC_MAX])
        plt.ylim([-SUCC_DISTANCE, SUCC_DISTANCE])

        ## Save figure
        plt.tight_layout()
        exportFigure(plt, plotOutDir, plotOutFileName, metadata)
        plt.close('all')

if __name__ == "__main__":
    Plot().cli()
