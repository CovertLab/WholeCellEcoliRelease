"""
Glucose mass yield distributions.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 9/19/19
"""

from __future__ import absolute_import, division, print_function

import cPickle
import os

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

from models.ecoli.analysis import variantAnalysisPlot
from models.ecoli.analysis.AnalysisPaths import AnalysisPaths
from models.ecoli.sim.variants.kinetic_constraints_factorial_experiments import get_disabled_constraints
from wholecell.analysis.analysis_tools import exportFigure
from wholecell.io.tablereader import TableReader
from wholecell.utils import constants, filepath, units
from wholecell.utils.sparkline import whitePadSparklineAxis


GLUCOSE_ID = 'GLC[p]'
FLUX_UNITS = units.mmol / units.g / units.h
MASS_UNITS = units.fg
GROWTH_UNITS = units.fg / units.s
ADDITIONAL_DISABLED_CONSTRAINTS = {
	'SUCCINATE-DEHYDROGENASE-UBIQUINONE-RXN-SUC/UBIQUINONE-8//FUM/CPD-9956.31.',
	'NADH-DEHYDROG-A-RXN-NADH/UBIQUINONE-8/PROTON//NAD/CPD-9956/PROTON.46. (reverse)',
	'INORGPYROPHOSPHAT-RXN[CCO-CYTOSOL]-PPI/WATER//Pi/PROTON.34.',
	'GLUTATHIONE-REDUCT-NADPH-RXN',
	}
VARIANT_LABELS = ['Succ/Fum\nDisabled', 'Succ/Fum\nEnabled', 'New Constraints']
N_VARIANTS = len(VARIANT_LABELS)
# Senior. Regulation of Nitrogen Metabolism in Escherichia coli and Klebsiella aerogenes:
# Studies with the Continuous-Culture Technique. 1975. Table 2 at growth rate of 0.9/hr (46 min doubling time)
VALIDATION_YIELD = 0.46

class Plot(variantAnalysisPlot.VariantAnalysisPlot):
	def do_plot(self, inputDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata):
		if not os.path.isdir(inputDir):
			raise Exception, 'inputDir does not currently exist as a directory'

		filepath.makedirs(plotOutDir)

		ap = AnalysisPaths(inputDir, variant_plot=True)
		all_variants = ap.get_variants()
		variants = -np.ones(N_VARIANTS)
		for v, variant in enumerate(all_variants):
			disable_constraints, additional_disabled = get_disabled_constraints(variant)
			if additional_disabled is None:
				variants[0] = variant
			elif len(additional_disabled) == 0:
				variants[1] = variant
			elif ADDITIONAL_DISABLED_CONSTRAINTS == set(additional_disabled):
				variants[2] = variant

		if np.any(variants < 0):
			print('Not enough variants to analyze')
			return

		with open(os.path.join(inputDir, 'kb', constants.SERIALIZED_FIT1_FILENAME), 'rb') as f:
			sim_data = cPickle.load(f)

		all_yields = []
		for variant in variants:
			yields = []

			for sim_dir in ap.get_cells(variant=[variant]):
				sim_out_dir = os.path.join(sim_dir, 'simOut')

				# Listeners used
				fba_reader = TableReader(os.path.join(sim_out_dir, 'FBAResults'))
				main_reader = TableReader(os.path.join(sim_out_dir, 'Main'))
				mass_reader = TableReader(os.path.join(sim_out_dir, 'Mass'))

				# Load data
				time_step_sec = main_reader.readColumn('timeStepSec')

				external_fluxes = fba_reader.readColumn('externalExchangeFluxes')
				external_molecules = fba_reader.readAttribute('externalMoleculeIDs')

				dry_mass = MASS_UNITS * mass_reader.readColumn('dryMass')
				growth = GROWTH_UNITS * mass_reader.readColumn('growth') / time_step_sec

				# Calculate growth yield on glucose
				glc_idx = external_molecules.index(GLUCOSE_ID)
				glc_flux = FLUX_UNITS * external_fluxes[:, glc_idx]
				glc_mw = sim_data.getter.getMass([GLUCOSE_ID])[0]
				glc_mass_flux = glc_flux * glc_mw * dry_mass
				glc_mass_yield = growth / -glc_mass_flux

				yields += list(glc_mass_yield[1:].asNumber())

			all_yields += [yields]

		for i, v1 in enumerate(variants):
			for j, v2 in enumerate(variants[i+1:]):
				t, p = stats.ttest_ind(all_yields[i], all_yields[i+j+1], equal_var=False)
				print('p={:.2e} for variant {} vs variant {}'.format(p, v1, v2))

		plt.figure(figsize=(4, 4))
		xticks = range(N_VARIANTS)

		# Plot data
		plt.violinplot(all_yields, xticks, showmeans=False, showextrema=False)
		plt.axhline(VALIDATION_YIELD, linestyle='--', color='#eb7037')

		# Format axes
		ax = plt.gca()
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		plt.xticks(xticks, VARIANT_LABELS)
		plt.ylabel('Glucose Yield\n(g cell / g glucose)')

		plt.tight_layout()
		exportFigure(plt, plotOutDir, plotOutFileName, metadata)
		plt.close('all')


if __name__ == '__main__':
	Plot().cli()
