'''
Analyze results from metabolism_kinetic_objective_weight variant
'''

from __future__ import absolute_import, division, print_function

import os
import re
from typing import Dict, List, Sequence, Tuple

from matplotlib import pyplot as plt
import numpy as np
from six.moves import cPickle, range

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import parallelization, units
from wholecell.utils.sparkline import whitePadSparklineAxis
from six.moves import zip


MODEL_FLUX_UNITS = COUNTS_UNITS / VOLUME_UNITS / TIME_UNITS
DCW_FLUX_UNITS = units.mmol / units.g / units.h

FRAC_CONC_OFF_AXIS = 0.05
FRAC_FLUX_OFF_AXIS = 0.05

OUTLIER_REACTIONS = [
	# Add reaction IDs to exclude from central carbon correlation
	]  # type: Sequence


def analyze_variant(args):
	# type: (Tuple[int, AnalysisPaths, List[str], np.ndarray, List[bool]]) -> tuple
	'''
	Function to analyze the data for each variant in parallel

	Inputs:
		variant (int) - variant number to analyze
		ap (AnalysisPaths object)
		toya_reactions (list of str) - toya reaction IDs
		toya_fluxes (array of floats) - flux associated with each reaction in toya_reactions
		outlier_filter (list of bool) - True if associated Toya reaction should be excluded
	'''

	variant, ap, toya_reactions, toya_fluxes, outlier_filter = args
	n_sims = 0

	# Load sim_data attributes for the given variant
	sim_data = cPickle.load(open(ap.get_variant_kb(variant), 'rb'))
	cell_density = sim_data.constants.cell_density
	n_avogadro = sim_data.constants.n_avogadro
	lambdas = sim_data.process.metabolism.kinetic_objective_weight

	# Lists for each cell in current variant
	growth_rate = []
	actual_conc = []
	target_conc = []
	homeostatic_objective_values = []
	kinetic_objective_values = []
	actual_flux = []
	target_flux = []
	toya_model_fluxes = {}  # type: Dict[str, List[np.ndarray]]
	for rxn in toya_reactions:
		toya_model_fluxes[rxn] = []

	for sim_dir in ap.get_cells(variant=[variant]):
		sim_out_dir = os.path.join(sim_dir, 'simOut')

		# Create readers for data
		try:
			mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))
			fba_results_reader = TableReader(os.path.join(sim_out_dir, 'FBAResults'))
			enzyme_kinetics_reader = TableReader(os.path.join(sim_out_dir, 'EnzymeKinetics'))
			bulk_reader = TableReader(os.path.join(sim_out_dir, 'BulkMolecules'))
		except Exception as e:
			print(e)
			continue

		# Mass related values
		try:
			cell_mass = units.fg * mass_reader.readColumn('cellMass')
			dry_mass = units.fg * mass_reader.readColumn('dryMass')
		except Exception as e:
			print(e)
			continue
		dcw_to_volume = cell_density * (dry_mass / cell_mass).asNumber()
		volume = cell_mass / cell_density

		# Growth rates
		# Growth rate stored in units of per second and first value will be nan
		growth = mass_reader.readColumn('instantaniousGrowthRate')
		if growth.size <= 1:
			continue
		growth_rate.append(np.mean(3600 * growth[1:]))

		# Metabolite comparison
		metabolite_ids = fba_results_reader.readAttribute('homeostaticTargetMolecules')
		bulk_ids = bulk_reader.readAttribute('objectNames')

		bulk_idxs = [bulk_ids.index(id_) for id_ in metabolite_ids]
		actual_counts = bulk_reader.readColumn('counts')[:, bulk_idxs]
		actual_conc.append(np.mean((1. / n_avogadro / volume * actual_counts.T).asNumber(COUNTS_UNITS / VOLUME_UNITS), axis=1))
		target_conc.append(np.nanmean(fba_results_reader.readColumn('targetConcentrations')[1:, :], axis=0))
		homeostatic_objective_values.append(np.mean(np.sum(fba_results_reader.readColumn('homeostaticObjectiveValues'), axis=1)))

		# Flux target comparison
		# Flux for first recorded step is 0
		target_fluxes = MODEL_FLUX_UNITS * enzyme_kinetics_reader.readColumn('targetFluxes').T
		actual_fluxes = MODEL_FLUX_UNITS * enzyme_kinetics_reader.readColumn('actualFluxes').T

		target_fluxes = (target_fluxes / dcw_to_volume).asNumber(DCW_FLUX_UNITS)
		actual_fluxes = (actual_fluxes / dcw_to_volume).asNumber(DCW_FLUX_UNITS)

		target_flux.append(np.nanmean(target_fluxes[:, 1:], axis=1))
		actual_flux.append(np.nanmean(actual_fluxes[:, 1:], axis=1))

		kinetic_objective = np.abs(1 - actual_fluxes / target_fluxes)
		filter_idx = ~np.isfinite(kinetic_objective)
		kinetic_objective[filter_idx] = 0
		kinetic_objective = np.sum(kinetic_objective, axis=0)
		kinetic_objective_values.append(np.mean(kinetic_objective))

		# Toya comparison
		# Toya units read in as mmol/g/hr
		reaction_ids = np.array(fba_results_reader.readAttribute('reactionIDs'))
		reaction_fluxes = MODEL_FLUX_UNITS * fba_results_reader.readColumn('reactionFluxes').T
		reaction_fluxes = (reaction_fluxes / dcw_to_volume).asNumber(DCW_FLUX_UNITS).T

		for toya_reaction_id in toya_reactions:
			flux_time_course = []  # type: List[np.ndarray]

			for rxn in reaction_ids:
				if re.findall(toya_reaction_id, rxn):
					reverse = 1
					if re.findall('(reverse)', rxn):
						reverse = -1

					if len(flux_time_course):
						flux_time_course += reverse * reaction_fluxes[1:, np.where(reaction_ids == rxn)]
					else:
						flux_time_course = reverse * reaction_fluxes[1:, np.where(reaction_ids == rxn)]

			if len(flux_time_course):
				flux_ave = np.mean(flux_time_course)
				toya_model_fluxes[toya_reaction_id].append(flux_ave)

		n_sims += 1

	# Growth rates
	growth_rates = np.mean(growth_rate)

	# Metabolite comparison
	actual_conc = np.mean(actual_conc, axis=0)
	target_conc = np.mean(target_conc, axis=0)
	conc_correlation = np.corrcoef(np.log(actual_conc), np.log(target_conc))[0, 1]
	n_conc_off_axis = np.sum(np.abs((target_conc - actual_conc) / target_conc) > FRAC_CONC_OFF_AXIS)

	# Flux target comparison
	# Add small value (1e-6) to fluxes for correlation so not taking log of 0
	actual_flux = np.mean(actual_flux, axis=0)
	target_flux = np.mean(target_flux, axis=0)
	flux_correlation = np.corrcoef(np.log(actual_flux + 1e-6), np.log(target_flux + 1e-6))[0, 1]
	n_flux_off_axis = np.sum(np.abs((target_flux - actual_flux) / target_flux) > FRAC_FLUX_OFF_AXIS)
	mask = (actual_flux != 0)
	nonzero_flux_correlation = np.corrcoef(np.log(actual_flux[mask] + 1e-6), np.log(target_flux[mask] + 1e-6))[0, 1]

	n_flux_above_0 = np.sum(actual_flux > 0) + np.sum((actual_flux == 0) & (target_flux == 0))

	# Toya comparison
	ave_toya_model = np.array([np.mean(toya_model_fluxes[rxn]) for rxn in toya_reactions])
	correlation_coefficient = np.corrcoef(ave_toya_model, toya_fluxes)[0, 1]
	filtered_correlation_coefficient = np.corrcoef(ave_toya_model[outlier_filter], toya_fluxes[outlier_filter])[0, 1]

	# Objective values
	# Need to filter nan and inf for kinetic
	kinetic_objective_value = np.mean(kinetic_objective_values)
	kinetic_objective_std = np.std(kinetic_objective_values)
	homeostatic_objective_value = np.mean(homeostatic_objective_values)
	homeostatic_objective_std = np.std(homeostatic_objective_values)

	n_metabolites = len(actual_conc)
	n_fluxes = len(actual_flux)

	return (lambdas, n_sims, growth_rates, conc_correlation, n_conc_off_axis, flux_correlation,
		n_flux_off_axis, nonzero_flux_correlation, n_flux_above_0, correlation_coefficient,
		filtered_correlation_coefficient, kinetic_objective_value, kinetic_objective_std,
		homeostatic_objective_value, homeostatic_objective_std, n_metabolites, n_fluxes)

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()
		n_variants = len(variants)
		total_sims = ap.n_seed * ap.n_generation

		if n_variants <= 1:
			print('This plot {} only runs for multiple variants'.format(__name__))
			return

		# Load validation data
		validation_data = cPickle.load(open(validationDataFile, 'rb'))
		toya_reactions = validation_data.reactionFlux.toya2010fluxes['reactionID']
		toya_fluxes = np.array([x.asNumber(DCW_FLUX_UNITS) for x in validation_data.reactionFlux.toya2010fluxes['reactionFlux']])
		outlier_filter = [False if rxn in OUTLIER_REACTIONS else True for rxn in toya_reactions]

		# Arrays to populate for plots
		lambdas = np.zeros(n_variants)
		n_sims = np.zeros(n_variants)
		growth_rates = np.zeros(n_variants)
		conc_correlation = np.zeros(n_variants)
		n_conc_off_axis = np.zeros(n_variants)
		flux_correlation = np.zeros(n_variants)
		nonzero_flux_correlation = np.zeros(n_variants)
		n_flux_above_0 = np.zeros(n_variants)
		n_flux_off_axis = np.zeros(n_variants)
		correlation_coefficient = np.zeros(n_variants)
		filtered_correlation_coefficient = np.zeros(n_variants)
		homeostatic_objective_value = np.zeros(n_variants)
		kinetic_objective_value = np.zeros(n_variants)
		homeostatic_objective_std = np.zeros(n_variants)
		kinetic_objective_std = np.zeros(n_variants)

		# Pull information from sim data and listeners in parallel
		pool = parallelization.pool(num_processes=self.cpus)
		args = list(zip(
			variants,
			[ap] * n_variants,
			[toya_reactions] * n_variants,
			[toya_fluxes] * n_variants,
			[outlier_filter] * n_variants
			))
		results = pool.map(analyze_variant, args)
		pool.close()
		pool.join()
		for i, result in enumerate(results):
			(lambdas[i],
				n_sims[i],
				growth_rates[i],
				conc_correlation[i],
				n_conc_off_axis[i],
				flux_correlation[i],
				n_flux_off_axis[i],
				nonzero_flux_correlation[i],
				n_flux_above_0[i],
				correlation_coefficient[i],
				filtered_correlation_coefficient[i],
				kinetic_objective_value[i],
				kinetic_objective_std[i],
				homeostatic_objective_value[i],
				homeostatic_objective_std[i],
				n_metabolites,
				n_fluxes) = result

		# Add each term of interest to an objective value to select the
		# lambda value that maximizes the objective
		scaled_growth_rate = growth_rates / growth_rates[0]
		scaled_growth_rate[scaled_growth_rate > 1] = 1
		objective = (
			n_sims / total_sims
			+ scaled_growth_rate
			+ conc_correlation
			+ (1 - n_conc_off_axis / n_metabolites)
			+ flux_correlation
			+ (1 - n_flux_off_axis / n_fluxes)
			+ nonzero_flux_correlation
			+ n_flux_above_0 / n_fluxes
			+ correlation_coefficient
		)

		tick_labels = [r'$10^{%i}$' % np.log10(x) if x != 0 else '0' for x in lambdas]
		lambdas = [np.log10(x) if x != 0 else np.nanmin(np.log10(lambdas[lambdas != 0]))-1 for x in lambdas]

		plt.figure(figsize = (8.5, 22))
		plt.style.use('seaborn-deep')
		subplots = 9

		# Growth rates
		ax = plt.subplot(subplots, 1, 1)
		plt.bar(lambdas, growth_rates / growth_rates[0], align='center')
		plt.axhline(1, linestyle='--', color='k')
		plt.ylim([0, 2])
		plt.ylabel('Growth rate deviation\nfrom no kinetics')
		whitePadSparklineAxis(ax, xAxis=False)
		plt.yticks([0, 1, 2])

		# Flux target comparisons
		ax = plt.subplot(subplots, 1, 2)
		plt.bar(lambdas, nonzero_flux_correlation, align='center')
		for lam, val in zip(lambdas, nonzero_flux_correlation):
			plt.text(lam, val, '{:.3f}'.format(val), ha='center')
		plt.ylim([0, 1])
		plt.ylabel('Kinetic target flux PCC')
		whitePadSparklineAxis(ax, xAxis=False)

		ax = plt.subplot(subplots, 1, 3)
		plt.bar(lambdas, n_flux_above_0 / n_fluxes, align='center')
		for lam, val in zip(lambdas, n_flux_above_0):
			plt.text(lam, val / n_fluxes, '{:.0f}/{:.0f}'.format(val, n_fluxes), ha='center')
		plt.ylim([0, 1])
		plt.ylabel('Fraction of fluxes\nabove 0')
		whitePadSparklineAxis(ax, xAxis=False)

		ax = plt.subplot(subplots, 1, 4)
		plt.bar(lambdas, n_flux_off_axis / n_fluxes, align='center')
		for lam, val in zip(lambdas, n_flux_off_axis):
			plt.text(lam, val / n_fluxes, '{:.0f}/{:.0f}'.format(val, n_fluxes), ha='center')
		plt.ylim([0, 1])
		plt.ylabel('Fraction of fluxes\noff axis (>{:.0f}%)'.format(FRAC_FLUX_OFF_AXIS*100))
		whitePadSparklineAxis(ax, xAxis=False)

		# Metabolite comparisons
		ax = plt.subplot(subplots, 1, 5)
		plt.bar(lambdas, conc_correlation, align='center')
		for lam, val in zip(lambdas, conc_correlation):
			plt.text(lam, val, '{:.3f}'.format(val), ha='center')
		plt.ylim([0, 1])
		plt.ylabel('Concentration PCC')
		whitePadSparklineAxis(ax, xAxis=False)

		ax = plt.subplot(subplots, 1, 6)
		plt.bar(lambdas, n_conc_off_axis / n_metabolites, align='center')
		for lam, val in zip(lambdas, n_conc_off_axis):
			plt.text(lam, val / n_metabolites, '{:.0f}/{:.0f}'.format(val, n_metabolites), ha='center')
		plt.ylim([0, 1])
		plt.ylabel('Fraction of concentrations\noff axis (>{:.0f}%)'.format(FRAC_CONC_OFF_AXIS*100))
		whitePadSparklineAxis(ax, xAxis=False)

		# Toya comparison
		ax = plt.subplot(subplots, 1, 7)
		plt.bar(lambdas, filtered_correlation_coefficient, align='center')
		for lam, val in zip(lambdas, filtered_correlation_coefficient):
			plt.text(lam, val, '{:.3f}'.format(val), ha='center')
		plt.ylim([0, 1])
		plt.ylabel('Central carbon flux PCC')
		whitePadSparklineAxis(ax, xAxis=False)

		# Viable sims
		ax = plt.subplot(subplots, 1, 8)
		plt.bar(lambdas, n_sims, align='center')
		for lam, val in zip(lambdas, n_sims):
			plt.text(lam, val, '{:.0f}'.format(val), ha='center')
		plt.ylabel('Number of sims\nwith data')
		whitePadSparklineAxis(ax, xAxis=False)

		# Lambda objective
		ax = plt.subplot(subplots, 1, 9)
		plt.bar(lambdas, objective, align='center')
		for lam, val in zip(lambdas, objective):
			plt.text(lam, val, '{:.3f}'.format(val), ha='center')
		plt.ylabel('Combined output objective')
		whitePadSparklineAxis(ax)
		plt.xticks(lambdas, tick_labels)

		plt.xlabel('lambda')

		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		# Plot kinetic vs homeostatic objective values
		plt.figure(figsize=(3.5, 3.5))
		ax = plt.gca()
		ax.set_xscale("log", nonposx='clip')
		ax.set_yscale("log", nonposy='clip')
		plt.errorbar(homeostatic_objective_value, kinetic_objective_value,
			xerr=homeostatic_objective_std, yerr=kinetic_objective_std,
			fmt='none', ecolor='k', alpha=0.5, linewidth=0.5)
		plt.plot(homeostatic_objective_value, kinetic_objective_value, "ob",
			markeredgewidth=0.1, alpha=0.9)
		for i in range(len(lambdas)):
			plt.text(homeostatic_objective_value[i], 0.6*kinetic_objective_value[i],
				i, horizontalalignment='center', verticalalignment='center')
		plt.xlabel('Homeostatic Objective Value')
		plt.ylabel('Kinetics Objective Value')

		whitePadSparklineAxis(ax)

		# Adjust limits to get tick labels to display
		xlim = ax.get_xlim()
		xlim = [10**np.floor(np.log10(xlim[0])), 10**np.ceil(np.log10(xlim[1]))]
		ax.set_xticks(xlim)
		ylim = ax.get_ylim()
		ylim = [10**np.floor(np.log10(ylim[0])), 10**np.ceil(np.log10(ylim[1]))]
		ax.set_yticks(ylim)

		exportFigure(plt, plotOutDir, '{}_obj'.format(plotOutFileName), metadata)

		plt.close('all')


if __name__ == "__main__":
	Plot().cli()
