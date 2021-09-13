"""
Analyzes parameters sensitivity from running variant param_sensitivity.
Outputs two plots showing sorted z score for each parameter's effect
on each output measure and individual parameter values for the most
significant parameters for each output difference measure.
"""

from __future__ import absolute_import, division, print_function

from functools import reduce
import io
import operator
import os
import re

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from scipy import special, stats
from six.moves import cPickle, range, zip
from typing import List, Optional, Tuple

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import COUNTS_UNITS, MASS_UNITS, TIME_UNITS, VOLUME_UNITS
from models.ecoli.sim.variants.param_sensitivity import number_params, split_indices
from reconstruction.ecoli.simulation_data import SimulationDataEcoli
from validation.ecoli.validation_data import ValidationDataEcoli
from wholecell.io import tsv
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import parallelization, sparkline, units


CONTROL_VARIANT = 0  # variant number for control simulation

ap = None  # type: Optional[AnalysisPaths]
sim_data = None  # type: Optional[SimulationDataEcoli]
validation_data = None  # type: Optional[ValidationDataEcoli]


def analyze_variant(args):
	# type: (Tuple[int, int]) -> np.ndarray
	'''
	Method to map each variant to for parallel analysis.

	Requires global variables:
		sim_data (SimulationData object)
		validation_data (ValidationData object)
		ap (AnalysisPaths object)

	Args:
		(variant, total_params): variant index;
			total number of parameters that are changed

	Returns:
		ndarray[float]: 2D array of results with each row corresponding to value below:
			number of times each parameter was increased
			number of times each parameter was decreased
			average growth rate for each parameter when increased
			average growth rate for each parameter when decreased
			average flux correlation for each parameter when increased
			average flux correlation for each parameter when decreased
	'''
	assert ap is not None
	assert sim_data is not None
	assert validation_data is not None

	variant, total_params = args
	if variant == 0:
		increase_indices = None
		decrease_indices = None
	else:
		increase_indices, decrease_indices = split_indices(sim_data, variant)

	increase_params_counts = np.zeros(total_params)
	decrease_params_counts = np.zeros(total_params)
	increase_params_growth_rate = np.zeros(total_params)
	decrease_params_growth_rate = np.zeros(total_params)
	increase_params_flux_correlation = np.zeros(total_params)
	decrease_params_flux_correlation = np.zeros(total_params)

	cell_density = sim_data.constants.cell_density
	flux_units = units.mmol / units.g / units.h

	# Validation data
	toya_reactions = validation_data.reactionFlux.toya2010fluxes["reactionID"]
	toya_fluxes = [f.asNumber(flux_units) for f in validation_data.reactionFlux.toya2010fluxes["reactionFlux"]]

	for sim_dir in ap.get_cells(variant=[variant]):
		simOutDir = os.path.join(sim_dir, "simOut")

		try:
			# Listeners used
			mass_reader = TableReader(os.path.join(simOutDir, 'Mass'))
			fba_results_reader = TableReader(os.path.join(simOutDir, "FBAResults"))

			# Load data
			## Growth rate
			growth_rate = np.nanmean(mass_reader.readColumn('instantaneous_growth_rate')[-5:]) * 3600  # 1/hr

			## Central carbon flux
			dry_mass = mass_reader.readColumn('dryMass')[-5:]
			cell_mass = mass_reader.readColumn('cellMass')[-5:]
			coefficient = dry_mass / cell_mass * cell_density.asNumber(MASS_UNITS / VOLUME_UNITS)

			reaction_ids = np.array(fba_results_reader.readAttribute("reactionIDs"))
			reaction_fluxes = (COUNTS_UNITS / MASS_UNITS / TIME_UNITS) * (fba_results_reader.readColumn("reactionFluxes")[-5:].T / coefficient).T
		# Exclude failed sims
		except Exception as e:
			print('Variant {} exception: {}'.format(variant, e))
			continue

		# Extract fluxes in Toya data set from simulation output
		model_fluxes = np.zeros_like(toya_fluxes)
		for i, toya_reaction in enumerate(toya_reactions):
			flux_time_course = []  # type: List[np.ndarray]

			for rxn in reaction_ids:
				if re.findall(toya_reaction, rxn):
					reverse = 1
					if re.findall("(reverse)", rxn):
						reverse = -1

					if len(flux_time_course):
						flux_time_course += reverse * reaction_fluxes[:, np.where(reaction_ids == rxn)]
					else:
						flux_time_course = reverse * reaction_fluxes[:, np.where(reaction_ids == rxn)]

			if len(flux_time_course):
				model_fluxes[i] = np.mean(flux_time_course).asNumber(flux_units)

		flux_r, _ = stats.pearsonr(toya_fluxes, model_fluxes)

		increase_params_counts[increase_indices] += 1
		decrease_params_counts[decrease_indices] += 1
		increase_params_growth_rate[increase_indices] += growth_rate
		decrease_params_growth_rate[decrease_indices] += growth_rate
		increase_params_flux_correlation[increase_indices] += flux_r
		decrease_params_flux_correlation[decrease_indices] += flux_r

	return np.vstack((increase_params_counts, decrease_params_counts,
		increase_params_growth_rate, decrease_params_growth_rate,
		increase_params_flux_correlation, decrease_params_flux_correlation,
		))

def headers(labels, name):
	'''
	Creates headers for tsv file

	Args:
		labels (list[str]): labels for each output value
		name (str): name of data type

	Returns:
		list[str]: combined names for header
	'''

	return ['{}\n{}'.format(name, label) for label in labels]


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if metadata.get('variant', '') != 'param_sensitivity':
			print('This plot only runs for the param_sensitivity variant.')
			return

		global ap
		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = np.array(ap.get_variants())

		# Check to analyze control (variant 0) separately from other variants
		use_control = False
		if CONTROL_VARIANT in variants:
			use_control = True
			variants = variants[variants != CONTROL_VARIANT]
		n_variants = len(variants)

		# Load one instance of sim_data to get number of parameters and ids
		global sim_data
		global validation_data
		with open(simDataFile, 'rb') as f:
			sim_data = cPickle.load(f)
		with open(validationDataFile, 'rb') as f:
			validation_data = cPickle.load(f)

		# sim_data information
		total_params = sum(number_params(sim_data))
		rna_to_gene = {gene['rna_ids'][0]: gene['symbol'] for gene in sim_data.process.replication.gene_data}
		monomer_to_gene = {gene['monomer_id']: gene['symbol'] for gene in sim_data.process.replication.gene_data}
		rna_ids = sim_data.process.transcription.rna_data['id']
		monomer_ids = sim_data.process.translation.monomer_data['id']

		# IDs must match order from param_indices() from param_sensitivity.py variant
		param_ids = np.array(
			['{} RNA deg rate'.format(rna_to_gene[rna[:-3]]) for rna in rna_ids]
			+ ['{} protein deg rate'.format(monomer_to_gene[monomer[:-3]]) for monomer in monomer_ids]
			+ ['{} translation eff'.format(monomer_to_gene[monomer[:-3]]) for monomer in monomer_ids]
			+ ['{} synth prob'.format(rna_to_gene[rna[:-3]]) for rna in rna_ids])
		if len(param_ids) != total_params:
			raise ValueError('Number of adjusted parameters and list of ids do not match.')

		pool = parallelization.pool()
		args = zip(
			variants,
			[total_params] * n_variants,
			)

		results = pool.imap_unordered(analyze_variant, args)
		(increase_params_counts,
			decrease_params_counts,
			increase_params_growth_rate,
			decrease_params_growth_rate,
			increase_params_flux_correlation,
			decrease_params_flux_correlation) = reduce(operator.add, results)
		pool.close()
		pool.join()

		# Calculate effects and z score
		labels = [
			'growth rate',
			'flux correlation',
			]
		increase_params_data = np.vstack((
			increase_params_growth_rate / increase_params_counts,
			increase_params_flux_correlation / increase_params_counts,
			))
		decrease_params_data = np.vstack((
			decrease_params_growth_rate / decrease_params_counts,
			decrease_params_flux_correlation / decrease_params_counts,
			))
		n_outputs = len(labels)

		# Difference between effect when parameter increased vs decreased
		data_diff = increase_params_data - decrease_params_data
		mean_diff = np.nanmean(data_diff, axis=1).reshape(-1, 1)
		std_diff = np.nanstd(data_diff, axis=1).reshape(-1, 1)
		z_score_diff = (data_diff - mean_diff) / std_diff

		# Individual increase or decrease effects to check asymmetric effects
		all_data = np.hstack((increase_params_data, decrease_params_data))
		mean = np.nanmean(all_data, axis=1).reshape(-1, 1)
		std = np.nanstd(all_data, axis=1).reshape(-1, 1)
		z_score_increase = (increase_params_data - mean) / std
		z_score_decrease = (decrease_params_data - mean) / std

		# Get control data
		if use_control:
			control_counts, _, control_growth_rate, _, control_flux_correlation, _ = analyze_variant((CONTROL_VARIANT, total_params))
			control_data = [
				control_growth_rate[0] / control_counts[0],
				control_flux_correlation[0] / control_counts[0],
				]
		else:
			control_data = [None] * n_outputs

		# Multiple hypothesis adjustment for significance of each parameter.
		# Solves Gaussian CDF for how many standard deviations are needed to
		# include 1 - 0.05 / total_params of the data (test each parameter for p<0.05).
		n_stds = special.erfinv(2 * (1 - 0.05 / total_params) - 1) * np.sqrt(2)

		# Plot histogram
		plt.figure(figsize=(10, 4*n_outputs))

		for i, (z_diff, z_increase, z_decrease) in enumerate(zip(z_score_diff, z_score_increase, z_score_decrease)):
			sorted_idx = np.argsort(z_diff)

			## Plot data
			ax = plt.subplot(n_outputs, 2, 2*i + 1)
			plt.yscale('symlog', linthreshold=0.01)
			plt.bar(list(range(total_params)), z_diff[sorted_idx])
			plt.axhline(n_stds , color='k', linestyle='--')
			plt.axhline(-n_stds, color='k', linestyle='--')

			## Format axes
			sparkline.whitePadSparklineAxis(ax, xAxis=False)
			plt.xticks([])
			plt.yticks([-n_stds, 0, n_stds])
			ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
			lim = np.max(np.abs(plt.ylim()))
			plt.ylim([-lim, lim])
			if i == 0:
				plt.title('Difference of Positive and Negative\nParameter Changes')
			if i == n_outputs - 1:
				plt.xlabel('Sorted Parameters')
			plt.ylabel('Z score\nparameter effect on {}\n(log scale)'.format(labels[i]))

			## Plot data
			ax = plt.subplot(n_outputs, 2, 2*i + 2)
			plt.yscale('symlog', linthreshold=0.01)
			plt.bar(list(range(total_params)), z_increase[sorted_idx], color='g')
			plt.bar(list(range(total_params)), z_decrease[sorted_idx], color='r')
			plt.axhline(n_stds , color='k', linestyle='--')
			plt.axhline(-n_stds, color='k', linestyle='--')

			## Format axes
			sparkline.whitePadSparklineAxis(ax, xAxis=False)
			plt.xticks([])
			plt.yticks([-n_stds, 0, n_stds])
			ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
			plt.ylim([-lim, lim])
			if i == 0:
				plt.title('Positive and Negative\nParameter Changes')
			if i == n_outputs - 1:
				plt.xlabel('Sorted Parameters')

		## Save figure
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)

		# Plot individual parameters
		individual_indices = [
			np.nanargmax(z_score_diff[0, :]),
			np.nanargmin(z_score_diff[0, :]),
			np.nanargmax(z_score_diff[1, :]),
			np.nanargmin(z_score_diff[1, :]),
			]
		n_individual = len(individual_indices)
		x_values = [-1, 0, 1]
		plt.figure()

		for i, label in enumerate(labels):
			shared_ax = None
			for j, idx in enumerate(individual_indices):
				## Shared y axis for each row
				ax = plt.subplot(n_outputs, n_individual, i*n_individual + j + 1, sharey=shared_ax)
				if shared_ax is None:
					shared_ax = ax

				## Plot data
				plt.plot(x_values, [decrease_params_data[i, idx], control_data[i], increase_params_data[i, idx]], 'x')

				## Format axes
				plt.xticks(x_values, ['Decrease', 'Control', 'Increase'])
				ax.tick_params(labelsize=6)
				ax.spines['right'].set_visible(False)
				ax.spines['top'].set_visible(False)
				if i < n_outputs - 1:
					ax.tick_params(labelbottom=False)
				if j > 0:
					ax.tick_params(labelleft=False)
				if i == 0:
					plt.title(param_ids[idx], fontsize=8)
				if j == 0:
					plt.ylabel(label, fontsize=7)

		## Save figure
		plt.tight_layout()
		exportFigure(plt, plotOutDir, '{}_individual'.format(plotOutFileName))
		plt.close('all')

		# Save z scores to tsv
		with io.open(os.path.join(plotOutDir, '{}.tsv'.format(plotOutFileName)), 'wb') as f:
			writer = tsv.writer(f)

			writer.writerow(
				['Parameter']
				+ headers(labels, 'Z-score, difference')
				+ headers(labels, 'Z-score, increase')
				+ headers(labels, 'Z-score, decrease')
				+ headers(labels, 'Raw average, difference')
				+ headers(labels, 'Raw average, increase')
				+ headers(labels, 'Raw average, decrease')
				)
			writer.writerows(np.hstack((
				param_ids.reshape(-1, 1),
				z_score_diff.T,
				z_score_increase.T,
				z_score_decrease.T,
				data_diff.T,
				increase_params_data.T,
				decrease_params_data.T
				)))


if __name__ == "__main__":
	Plot().cli()
