# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"KmOptimization.py",
	"KmRNAdecayComparison.py",
	"aaCounts.py",
	"aaExchangeFluxes.py",
	"active_rnap_coordinates.py",
	"allReactionFluxes.py",
	"centralCarbonMetabolism.py",
	"centralCarbonMetabolismCorrelationTimeCourse.py",
	"centralCarbonMetabolismScatter.py",
	"concentrationDeviation.py",
	"dnaA_box_dynamics.py",
	"dntpCounts.py",
	"equilibriumComparison.py",
	"evaluationTime.py",
	"expression_rna_01_low.py",
	"expression_rna_02_med.py",
	"expression_rna_03_high.py",
	"fbaOptimizationProblem.py",
	"flagella_expression.py",
	"glucoseAndOxygenExchangeFluxes.py",
	"glucoseMassYield.py",
	"growthLimits.py",
	"inter_rnap_distance.py",
	"kineticsFluxComparison.py",
	"kineticsFluxComparisonKcatOnly.py",
	"mRnaHalfLives.py",
	"massFractionSummary.py",
	"massFractions.py",
	"metaboliteComparison.py",
	"metabolites.py",
	"mrnaCounts.py",
	"mrnaVsProteinCounts.py",
	"ntpCounts.py",
	"processMassBalance.py",
	"processMassBalanceDynamics.py",
	"proteinCounts.py",
	"proteinCountsValidation.py",
	"replication.py",
	"replisome_counts.py",
	"replisome_rnap_collisions.py",
	"ribosome30SCounts.py",
	"ribosome50SCounts.py",
	"ribosomeCapacity.py",
	"ribosomeCounts.py",
	"rnaDegradationCounts.py",
	"rnaSynthesisProbabilities.py",
	"rnapCapacity.py",
	"rnapCounts.py",
	"rnaseCounts.py",
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
		"dntpCounts.py",
		"evaluationTime.py",
		"glucoseAndOxygenExchangeFluxes.py",
		"massFractionSummary.py",
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
		"trna_charging.py",
		],
	'METABOLISM': [
		"aaExchangeFluxes.py",
		"allReactionFluxes.py",
		"centralCarbonMetabolism.py",
		"centralCarbonMetabolismCorrelationTimeCourse.py",
		"centralCarbonMetabolismScatter.py",
		"concentrationDeviation.py",
		"fbaOptimizationProblem.py",
		"glucoseAndOxygenExchangeFluxes.py",
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
		"transient_gene_dosage.py",
		"trpRegulation.py",
		"twoComponentSystem.py",
		],
	'REPLICATION': [
		"dnaA_box_dynamics.py",
		"replication.py",
		"replisome_counts.py",
		"transient_gene_dosage.py",
		"active_rnap_coordinates.py",
		"replisome_rnap_collisions.py",
		],
	'TRANSCRIPTION': [
		"inter_rnap_distance.py",
		"mRnaHalfLives.py",
		"mrnaCounts.py",
		"mrnaVsProteinCounts.py",
		"rnaSynthesisProbabilities.py",
		"rnapCapacity.py",
		"rnapCounts.py",
		"tRnaCounts.py",
		],
	'TRANSLATION': [
		"growthLimits.py",
		"mrnaVsProteinCounts.py",
		"proteinCounts.py",
		"ribosome30SCounts.py",
		"ribosome50SCounts.py",
		"ribosomeCapacity.py",
		"ribosomeCounts.py",
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
