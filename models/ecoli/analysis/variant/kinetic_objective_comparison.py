"""
Violin plots to compare predicted k_cats of the new list of disabled constraints with the baseline list.

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS
from models.ecoli.sim.variants.kinetic_constraints_factorial_experiments import get_disabled_constraints
from wholecell.analysis.analysis_tools import exportFigure, read_bulk_molecule_counts
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, filepath, units


# additional disabled constraints that are to be compared to baseline
ADDITIONAL_DISABLED_CONSTRAINTS = {
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)',
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.',
	'GLUTATHIONE-REDUCT-NADPH-RXN',
	}

OLD_MEASUREMENTS = {
	'ISOCITDEH-RXN': {
		'measurements': [106.3, 88.1], 'temps': [40, 40]},
	'GLYOXYLATE-REDUCTASE-NADP+-RXN__CPLX0-235': {
		'measurements': [203], 'temps': [25]},
	'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.': {
		'measurements': [1128, 177, 250, 230], 'temps': [38, 30, 30, 30]},
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.': {
		'measurements': [78, 110, 24, 85], 'temps': [30, 30, 30, 30]},
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)': {
		'measurements': [26], 'temps': [30]},
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.': {
		'measurements': [187, 390, 390], 'temps': [25, 25, 25]},
	'CYTDEAM-RXN': {
		'measurements': [185, 165, 49.68 ,132 ,45], 'temps': [25, 37, 25, 30, 25]},
	'GLUTATHIONE-REDUCT-NADPH-RXN': {
		'measurements': [733.3], 'temps': [30]},
	'PSERTRANSAM-RXN': {
		'measurements': [1.75], 'temps': [37]},
	'CITSYN-RXN__CITRATE-SI-SYNTHASE': {
		'measurements': [81], 'temps': [25]},
	}

NEW_MEASUREMENTS = {
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)': {
		'measurements': [600], 'temps': [30]},
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.': {
		'measurements': [42], 'temps': [25]},
	}

# Taken from paper/kinetic_constraints.tsv
# KM in units of uM
SIMULATION_KMS = {
	'ISOCITDEH-RXN': {
		'metabolite': 'NADP[c]', 'KM': 39.2, 'constraint_index': 188},
	'R601-RXN-FUM/REDUCED-MENAQUINONE//SUC/CPD-9728.38.': {
		'metabolite': 'FUM[c]', 'KM': 4, 'constraint_index': 262},
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.': {
		'metabolite': 'SUC[c]', 'KM': 2, 'constraint_index': 394},
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)': {
		'metabolite': 'NADH[c]', 'KM': 13, 'constraint_index': 222},
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.': {
		'metabolite': 'PPI[c]', 'KM': 0.13, 'constraint_index': 183},
	'CYTDEAM-RXN': {
		'metabolite': 'CYTOSINE[c]', 'KM': 200, 'constraint_index': 87},
	'CITSYN-RXN__CITRATE-SI-SYNTHASE': {
		'metabolite': 'ACETYL-COA[c]', 'KM': 120, 'constraint_index': 81},
	}

REACTIONS = sorted(OLD_MEASUREMENTS.keys())

# ignore data from first time steps
START_TIME_STEP = 2

def set_ticks(ax, labels):
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().set_tick_params(direction='out')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks(np.arange(1, len(labels) + 1))
	ax.set_xticklabels(labels, fontsize=8)
	ax.set_xlim(0.25, len(labels) + 0.75)


class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		variants = ap.get_variants()

		# scan all variants to find variant indexes for comparison
		old_variant = None
		new_variant = None
		for v, variant in enumerate(variants):
			disable_constraints, additional_disabled = get_disabled_constraints(variant)
			if additional_disabled is None:
				old_variant = variant
			elif ADDITIONAL_DISABLED_CONSTRAINTS == set(additional_disabled):
				new_variant = variant

		# if the baseline variant or the new variant are missing, stop plotting
		if (old_variant is None) or (new_variant is None):
			print('Variant simulations missing!')
			return

		compared_variants = [old_variant, new_variant]

		# Load sim_data
		with open(os.path.join(inputDir, 'kb', constants.SERIALIZED_FIT1_FILENAME), 'rb') as f:
			sim_data = cPickle.load(f)

		# get reactions from sim_data
		reactionCatalysts = sim_data.process.metabolism.reactionCatalysts

		reaction_to_enzyme = {r: reactionCatalysts[r][0] for r in REACTIONS}
		enzyme_names = reaction_to_enzyme.values()
		reactions_with_km = sorted(SIMULATION_KMS)
		km_metabolites = [SIMULATION_KMS[r]['metabolite'] for r in reactions_with_km]
		kms = np.array([SIMULATION_KMS[r]['KM'] for r in reactions_with_km])
		km_constraint_indices = [SIMULATION_KMS[r]['constraint_index'] for r in reactions_with_km]

		# initialize dictionaries for fluxes and concentrations
		all_reaction_fluxes = {}
		all_enzyme_concentrations = {}
		all_km_adjustments = {}
		for variant in compared_variants:
			reaction_fluxes = {r: [] for r in REACTIONS}
			enzyme_concentrations = {e: [] for e in enzyme_names}
			km_adjustments = {r: [] for r in reactions_with_km}
			for sim_dir in ap.get_cells(variant=[variant]):
				simOutDir = os.path.join(sim_dir, "simOut")

				# Listeners used
				try:
					kinetics_reader = TableReader(os.path.join(simOutDir, 'EnzymeKinetics'))
					fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
				except Exception as e:
					print(e)
					continue

				# read from kinetics listener
				counts_to_molar = ((COUNTS_UNITS / VOLUME_UNITS)
					* kinetics_reader.readColumn('countsToMolar')[START_TIME_STEP:].reshape(-1, 1)
					)
				all_constraints_used = kinetics_reader.readColumn('reactionConstraint')[START_TIME_STEP:]

				# Store fluxes
				reactionIDs = np.array(fbaResults.readAttribute("reactionIDs"))
				reactionFluxes = fbaResults.readColumn("reactionFluxes")[START_TIME_STEP:, :]
				reaction_flux_dict = dict(zip(reactionIDs, reactionFluxes.T))
				for reaction_id in REACTIONS:
					reaction_fluxes[reaction_id].extend(list(reaction_flux_dict[reaction_id]))

				# Store enzyme concentrations
				enzyme_counts, met_counts = read_bulk_molecule_counts(simOutDir, (enzyme_names, km_metabolites))
				enzyme_conc = counts_to_molar.asNumber(COUNTS_UNITS / VOLUME_UNITS) * enzyme_counts[START_TIME_STEP:, :]
				met_conc = counts_to_molar.asNumber(units.umol / units.L) * met_counts[START_TIME_STEP:, :]
				for enzyme_id, conc_time_series in zip(enzyme_names, enzyme_conc.T):
					enzyme_concentrations[enzyme_id].extend(list(conc_time_series))

				# Calculate enzyme saturation for reactions with KM values
				adjust_km = np.zeros((len(counts_to_molar), len(km_constraint_indices)), bool)
				for i, idx in enumerate(km_constraint_indices):
					constraint_used, _ = np.where(all_constraints_used == idx)
					adjust_km[constraint_used, i] = True
				enzyme_saturation = met_conc / (met_conc + kms)
				enzyme_saturation[~adjust_km] = 1
				for rxn, saturation in zip(reactions_with_km, enzyme_saturation.T):
					km_adjustments[rxn].extend(list(saturation))

			all_reaction_fluxes[variant] = reaction_fluxes
			all_enzyme_concentrations[variant] = enzyme_concentrations
			all_km_adjustments[variant] = km_adjustments

		### Make figure ###
		cols = 1
		rows = len(REACTIONS)
		plt.figure(figsize=(cols * 3, rows * 5))

		# go through each reaction to show predicted k_cat distribution for the
		# new and old variant, and experimental measurements
		for reaction_idx, reaction_id in enumerate(REACTIONS):
			enzyme_id = reaction_to_enzyme[reaction_id]

			# old measurements
			reaction_measurements = OLD_MEASUREMENTS[reaction_id]
			measurements = reaction_measurements['measurements']
			temps = reaction_measurements['temps']
			adjusted_measurements = np.array([2**((37.-t)/10.)*m for (m, t) in zip(measurements, temps)])

			# new measurements
			reaction_measurements = NEW_MEASUREMENTS.get(reaction_id, {})
			measurements = reaction_measurements.get('measurements', [])
			temps = reaction_measurements.get('temps', [])
			new_adjusted_measurements = np.array([2**((37.-t)/10.)*m for (m, t) in zip(measurements, temps)])

			# get effective kcat for GLUTATHIONE-REDUCT
			if reaction_id == 'GLUTATHIONE-REDUCT-NADPH-RXN':
				# saturated_fraction calculated from Smirnova, et al. (2005). "Effects of cystine and
				# hydrogen peroxideon glutathione status and expression of	antioxidant	genes in Escherichia coli"
				# Oxidized glutathione (GSSG in table 2) gives ~19 uM concentration (with 0.3 dry fraction and 1.1 g/mL density)
				# With 61 uM Km for this reaction, that gives a saturated fraction of 0.238
				saturated_fraction = 0.238
				new_adjusted_measurements = adjusted_measurements * saturated_fraction

			# Initialize subplots
			ax = plt.subplot(rows, cols, reaction_idx+1)

			# calculate the reaction's k_cat distribution for each compared variant
			k_cat_distribution = {}
			for variant in compared_variants:
				## Get data
				rxn_fluxes = np.array(all_reaction_fluxes[variant][reaction_id]) 			# mmol / L / s
				enzyme_concs = np.array(all_enzyme_concentrations[variant][enzyme_id])  	# mmol / L
				saturation = np.array(all_km_adjustments[variant].get(reaction_id, [1] * len(rxn_fluxes)))

				# calculate k_cats (adjusted for saturation in the sim), remove zeros, save to this variant's distribution
				k_cats = rxn_fluxes / enzyme_concs / saturation
				k_cats = k_cats[k_cats > 1e-10]
				k_cat_distribution[variant] = k_cats

			data = [k_cat_distribution[old_variant], k_cat_distribution[new_variant]]

			# plot
			violin_pos = [1, 3]	# position of violin plots [old, new]
			measure_pos = 2  	# position of measurements
			ax.violinplot(data, violin_pos, widths=1.0, showmeans=False, showextrema=False, showmedians=False)
			ax.scatter(np.full_like(adjusted_measurements, measure_pos), adjusted_measurements, marker='o', color='#eb7037', s=50, alpha=0.7)
			ax.scatter(np.full_like(new_adjusted_measurements, measure_pos), new_adjusted_measurements, marker='o', color='#eb7037', s=50, alpha=0.7)

			# format
			rxn_id_length = 25
			text_reaction_id = ('reaction: %s' % reaction_id[:rxn_id_length])
			labels = ['\nModel Predicted\n(Old Constraints)', 'Measured', '\nModel Predicted\n(New Constraints)']
			ax.set_title(text_reaction_id, fontsize=8)
			ax.set_ylabel('$k_{cat}$ (1/s)', fontsize=8)
			set_ticks(ax, labels)
			ax.set_yscale('log')

		### Create Plot ###
		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == "__main__":
	Plot().cli()