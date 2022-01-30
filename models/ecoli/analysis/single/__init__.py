from __future__ import absolute_import, division, print_function

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"KmOptimization.py",
	"KmRNAdecayComparison.py",
	"aaCounts.py",
	"aaExchangeFluxes.py",
	"active_rnap_coordinates.py",
	"all_enzyme_absences.py",
	"cell_wall_expression.py",
	"centralCarbonMetabolism.py",
	"centralCarbonMetabolismCorrelationTimeCourse.py",
	"centralCarbonMetabolismScatter.py",
	"chromosome_visualization.py",
	"concentrationDeviation.py",
	"compartment_mass_fraction_summary.py",
	"cotranscriptional_translation.py",
	"dnaA_box_dynamics.py",
	"dntpCounts.py",
	"equilibriumComparison.py",
	"evaluationTime.py",
	"expected_mechanistic_vs_real_uptake.py",
	"expression_rna_01_low.py",
	"expression_rna_02_med.py",
	"expression_rna_03_high.py",
	"external_exchange_fluxes.py",
	"fbaOptimizationProblem.py",
	"flagella_expression.py",
	"flux_bounds.py",
	"glucoseMassYield.py",
	"growthLimits.py",
	"inter_rnap_distance.py",
	"kineticsFluxComparison.py",
	"kineticsFluxComparisonKcatOnly.py",
	"mRnaHalfLives.py",
	"massFractionSummary.py",
	"mass_fractions.py",
	"mass_fractions_voronoi.py",
	"metaboliteComparison.py",
	"metabolites.py",
	"mrnaCounts.py",
	"mrnaVsProteinCounts.py",
	"ntpCounts.py",
	"outer_membrane_expression.py",
	"processMassBalance.py",
	"processMassBalanceDynamics.py",
	"proteinCounts.py",
	"proteinCountsValidation.py",
	"replication.py",
	"replisome_counts.py",
	"replisome_rnap_collisions.py",
	"replisome_rnap_collisions_by_gene.py",
	"ribosome_limitation.py",
	"ribosome30SCounts.py",
	"ribosome50SCounts.py",
	"ribosomeCapacity.py",
	"ribosomeCounts.py",
	"rnaDegradationCounts.py",
	"rnaSynthesisProbabilities.py",
	"rnapCapacity.py",
	"rnapCounts.py",
	"rnap_stalled.py",
	"rnaseCounts.py",
	"superhelical_density.py",
	"surface_area_comparison.py",
	"transcriptional_attenuation.py",
	"transient_gene_dosage.py",
	"trna_charging.py",
	"tRnaCounts.py",
	"trpRegulation.py",
	"twoComponentSystem.py",
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"aaCounts.py",
		"compartment_mass_fraction_summary.py",
		"dntpCounts.py",
		"evaluationTime.py",
		"external_exchange_fluxes.py",
		"massFractionSummary.py",
		"metabolites.py",
		"mrnaCounts.py",
		"ntpCounts.py",
		"processMassBalance.py",
		"processMassBalanceDynamics.py",
		"proteinCounts.py",
		"replication.py",
		"ribosomeCapacity.py",
		"ribosomeCounts.py",
		"rnapCapacity.py",
		"rnapCounts.py",
		"rnaseCounts.py",
		],
	'PARCA': [
		"KmOptimization.py",
		"KmRNAdecayComparison.py",
		],
	'GROWTH': [
		"growthLimits.py",  # TODO(jerry): include this?
		"ribosome_limitation.py",
		"trna_charging.py",
		],
	'KINETICS': [
		"centralCarbonMetabolismScatter.py",
		"expected_mechanistic_vs_real_uptake.py",
		"external_exchange_fluxes.py",
		"kineticsFluxComparison.py",
		],
	'METABOLISM': [
		"aaExchangeFluxes.py",
		"all_enzyme_absences.py",
		"centralCarbonMetabolism.py",
		"centralCarbonMetabolismCorrelationTimeCourse.py",
		"centralCarbonMetabolismScatter.py",
		"concentrationDeviation.py",
		"expected_mechanistic_vs_real_uptake.py",
		"external_exchange_fluxes.py",
		"fbaOptimizationProblem.py",
		"flux_bounds.py",
		"glucoseMassYield.py",
		"kineticsFluxComparison.py",
		"kineticsFluxComparisonKcatOnly.py",
		"metaboliteComparison.py",
		"metabolites.py",
		],
	'PAPER': [
		"expression_rna_01_low.py",
		"expression_rna_02_med.py",
		"expression_rna_03_high.py",
		],
	'RNA_DEGRADATION': [
		"rnaDegradationCounts.py",
		"rnaseCounts.py",
		],
	'REGULATION': [
		"equilibriumComparison.py",
		"transcriptional_attenuation.py",
		"transient_gene_dosage.py",
		"trpRegulation.py",
		"twoComponentSystem.py",
		],
	'REPLICATION': [
		"active_rnap_coordinates.py",
		"chromosome_visualization.py",
		"dnaA_box_dynamics.py",
		"replication.py",
		"replisome_counts.py",
		"transient_gene_dosage.py",
		"replisome_rnap_collisions.py",
		"replisome_rnap_collisions_by_gene.py",
		],
	'TRANSCRIPTION': [
		"cotranscriptional_translation.py",
		"inter_rnap_distance.py",
		"mRnaHalfLives.py",
		"mrnaCounts.py",
		"mrnaVsProteinCounts.py",
		"rnaSynthesisProbabilities.py",
		"rnapCapacity.py",
		"rnapCounts.py",
		"superhelical_density.py",
		"tRnaCounts.py",
		],
	'TRANSLATION': [
		"cotranscriptional_translation.py",
		"growthLimits.py",
		"mrnaVsProteinCounts.py",
		"proteinCounts.py",
		"ribosome_limitation.py",
		"ribosome30SCounts.py",
		"ribosome50SCounts.py",
		"ribosomeCapacity.py",
		"ribosomeCounts.py",
		"transcriptional_attenuation.py",
		"trna_charging.py",
		],
	'TWO_COMPONENT_SYSTEM': [
		"twoComponentSystem.py",
		],
	'VALIDATION': [
		"centralCarbonMetabolism.py",
		"centralCarbonMetabolismCorrelationTimeCourse.py",
		"centralCarbonMetabolismScatter.py",
		"inter_rnap_distance.py",
		"proteinCountsValidation.py",
		],
	}
