# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"biosynthesisProductionDynamics.py",
	"carRegulation.py",
	"cellCycleLength.py",
	"centralCarbonMetabolismCorrelationTimeCourse.py",
	"centralCarbonMetabolismScatter.py",
	"charging_molecules.py",
	"environmental_shift_fluxes.py",
	"functionalUnits.py",
	"functionalUnitsFC.py",
	"growthAffectingPolymerases.py",
	# "growthRateControl.py",
	# "initiationNoiseVsFoldChange.py",
	"kineticsFluxComparison.py",
	"limitedEnzymeGlutcyslig.py",
	"limitedEnzymeUgd.py",
	"limitedMetabolites.py",
	"massFractionSummary.py",
	"massFractionToUnity.py",
	"massShift.py",
	"mene_limitations.py",
	# "probProteinExistAndDouble.py",
	"proteinAvgCountVsBurstSize.py",
	"proteinConcentrations.py",
	"proteinCountVsFoldChange.py",
	"proteinCountsValidation.py",
	"proteinExistVsBurstSize.py",
	"proteinFoldChangeVsRnaDeg.py",
	"proteinFoldChangeVsTranscriptionEvents.py",
	"proteinFoldChangeVsTranslationEff.py",
	"proteomeMassFractions.py",
	"replication.py",
	"ribosomeProduction.py",
	"ribosomeUsage.py",
	# "rnaVsProteinPerCell.py",
	"rna_decay_03_high.py",
	"subgenerationalTranscription.py",
	"timeStep.py",
	"transcriptionEvents.py",
	"transcriptionFrequency.py",
	"transcriptionFrequencyOrdered.py",
	# "transcriptionGenomeCoverage.py",
	"translationFrequency.py",
	"trpRegulation.py",
	"trpSynthaseCapacityVsUsage.py",
	"tyrRegulation.py",
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"massFractionSummary.py",
		"massFractionToUnity.py",
		"replication.py",
		"timeStep.py",
		],
	'DIVISION': [
		"cellCycleLength.py",
		],
	'ENVIRONMENT_SHIFT': [
		"environmental_shift_fluxes.py",  # TODO(jerry): include this?
		"massShift.py",
		],
	'GROWTH': [
		"charging_molecules.py",
		"growthAffectingPolymerases.py",
		# "growthRateControl.py",
		],
	'METABOLISM': [
		"centralCarbonMetabolismCorrelationTimeCourse.py",
		"centralCarbonMetabolismScatter.py",
		"kineticsFluxComparison.py",
		"limitedMetabolites.py",
		],
	'PAPER': [
		"environmental_shift_fluxes.py",
		"functionalUnits.py",
		"functionalUnitsFC.py",
		"limitedEnzymeGlutcyslig.py",
		"limitedEnzymeUgd.py",
		"mene_limitations.py",
		"proteinAvgCountVsBurstSize.py",
		"proteinCountVsFoldChange.py",
		"proteinCountsValidation.py",
		"proteinExistVsBurstSize.py",
		"proteinFoldChangeVsRnaDeg.py",
		"transcriptionFrequencyOrdered.py",
		"trpRegulation.py",
		],
	'RNA_DEGRADATION': [
		"rna_decay_03_high.py",
		],
	'REGULATION': [
		"biosynthesisProductionDynamics.py",
		"carRegulation.py",
		"trpRegulation.py",
		"trpSynthaseCapacityVsUsage.py",
		"tyrRegulation.py",
		],
	'SUBGENERATIONAL': [
		"proteinAvgCountVsBurstSize.py",
		"proteinCountVsFoldChange.py",
		"proteinExistVsBurstSize.py",
		"proteinFoldChangeVsRnaDeg.py",
		"proteinFoldChangeVsTranscriptionEvents.py",
		"proteinFoldChangeVsTranslationEff.py",
		"subgenerationalTranscription.py",
		"transcriptionEvents.py",
		"transcriptionFrequency.py",
		"transcriptionFrequencyOrdered.py",
		"translationFrequency.py",
		],
	'TRANSLATION': [
		"charging_molecules.py",
		"proteinFoldChangeVsTranslationEff.py",  # TODO(jerry): include this?
		"ribosomeProduction.py",
		"ribosomeUsage.py",
		"translationFrequency.py",  # TODO(jerry): include this?
		],
	'VALIDATION': [
		"proteinCountsValidation.py",
		]
	}
