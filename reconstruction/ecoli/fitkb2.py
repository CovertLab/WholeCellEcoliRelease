
from __future__ import division
import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.compendium import growth_data

from wholecell.utils import units
from wholecell.utils.fitting import countsFromMassAndExpression, calcProteinCounts, calcProteinDistribution, calcProteinTotalCounts

from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

N_SEEDS = 20

def fitKb_2(kb, simOutDir):

	growthData = growth_data.GrowthData(kb)
	massFractions60 = growthData.massFractions(60)
	proteinMass = massFractions60["proteinMass"].asUnit(units.g)
	rnaMass = massFractions60["rnaMass"].asUnit(units.g)

	# Construct bulk container

	# We want to know something about the distribution of the copy numbers of
	# macromolecules in the cell.  While RNA and protein expression can be
	# approximated using well-described statistical distributions, we need
	# absolute copy numbers to form complexes.  To get a distribution, we must
	# instantiate many cells, form complexes, and finally compute the
	# statistics we will use in the fitting operations.

	bulkContainer = BulkObjectsContainer(kb.state.bulkMolecules.bulkData['id'])
	rnaView = bulkContainer.countsView(kb.process.transcription.rnaData["id"])
	proteinView = bulkContainer.countsView(kb.process.translation.monomerData["id"])
	complexationMoleculesView = bulkContainer.countsView(kb.process.complexation.moleculeNames)
	allMoleculesIDs = list(
		set(kb.process.transcription.rnaData["id"]) | set(kb.process.translation.monomerData["id"]) | set(kb.process.complexation.moleculeNames)
		)
	allMoleculesView = bulkContainer.countsView(allMoleculesIDs)

	allMoleculeCounts = np.empty((N_SEEDS, allMoleculesView.counts().size), np.int64)

	complexationStoichMatrix = kb.process.complexation.stoichMatrix().astype(np.int64, order = "F")

	complexationPrebuiltMatrices = mccBuildMatrices(
		complexationStoichMatrix
		)

	rnaDistribution = kb.process.transcription.rnaData['expression']

	rnaTotalCounts = countsFromMassAndExpression(
		rnaMass.asNumber(units.g),
		kb.process.transcription.rnaData["mw"].asNumber(units.g / units.mol),
		rnaDistribution,
		kb.constants.nAvogadro.asNumber(1 / units.mol)
		)

	proteinDistribution = calcProteinDistribution(kb)

	proteinTotalCounts = calcProteinTotalCounts(kb, proteinMass, proteinDistribution)

	for seed in xrange(N_SEEDS):
		randomState = np.random.RandomState(seed)

		allMoleculesView.countsIs(0)

		rnaView.countsIs(randomState.multinomial(
			rnaTotalCounts,
			rnaDistribution
			))

		proteinView.countsIs(randomState.multinomial(
			proteinTotalCounts,
			proteinDistribution
			))

		complexationMoleculeCounts = complexationMoleculesView.counts()

		updatedCompMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
			complexationMoleculeCounts,
			seed,
			complexationStoichMatrix,
			*complexationPrebuiltMatrices
			)

		complexationMoleculesView.countsIs(updatedCompMoleculeCounts)

		allMoleculeCounts[seed, :] = allMoleculesView.counts()

	bulkAverageContainer = BulkObjectsContainer(kb.state.bulkMolecules.bulkData['id'], np.float64)
	bulkDeviationContainer = BulkObjectsContainer(kb.state.bulkMolecules.bulkData['id'], np.float64)

	bulkAverageContainer.countsIs(allMoleculeCounts.mean(0), allMoleculesIDs)
	bulkDeviationContainer.countsIs(allMoleculeCounts.std(0), allMoleculesIDs)

	# Free up memory
	# TODO: make this more functional; one function for returning average & distribution
	del allMoleculeCounts
	del bulkContainer

	# ----- tRNA synthetase turnover rates ------
	# Fit tRNA synthetase kcat values based on expected rates of translation
	# compute values at initial time point

	## Compute rate of AA incorperation
	proteinComposition = kb.process.translation.monomerData["aaCounts"]
	initialDryMass = kb.constants.avgCellDryMassInit

	proteinMassFraction = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].asNumber(units.min) == 60.0
		]["proteinMassFraction"]

	initialProteinMass = initialDryMass * proteinMassFraction

	initialProteinCounts = calcProteinCounts(kb, initialProteinMass)

	initialProteinTranslationRate = (
		(np.log(2) / kb.constants.cellCycleLen + kb.process.translation.monomerData["degRate"]) * initialProteinCounts
		).asUnit(1 / units.s)

	initialAAPolymerizationRate = units.dot(
		units.transpose(proteinComposition), initialProteinTranslationRate
		).asUnit(units.aa / units.s)

	## Compute expression of tRNA synthetases
	## Assuming independence in variance
	synthetase_counts_by_group = np.zeros(len(kb.aa_synthetase_groups), dtype = np.float64)
	synthetase_variance_by_group = np.zeros(len(kb.aa_synthetase_groups), dtype = np.float)
	for idx, synthetase_group in enumerate(kb.aa_synthetase_groups.itervalues()):
		group_count = 0.
		group_variance = 0.
		for synthetase in synthetase_group:
			counts = bulkAverageContainer.countsView([synthetase]).counts()
			variance = bulkDeviationContainer.countsView([synthetase]).counts()
			group_count += counts
			group_variance += variance
		synthetase_counts_by_group[idx] = group_count
		synthetase_variance_by_group[idx] = group_variance

	## Saved for plotting
	kb.synthetase_counts = synthetase_counts_by_group
	kb.synthetase_variance = synthetase_variance_by_group
	kb.initial_aa_polymerization_rate = initialAAPolymerizationRate
	kb.minimum_trna_synthetase_rates = initialAAPolymerizationRate / synthetase_counts_by_group

	# TODO: Reimplement this with better fit taking into account the variance in aa
	#		utilization.
	## Scaling synthetase counts by -2*variance so that rates will be high enough
	## to accomodate stochastic behavior in the model without translation stalling.
	# scaled_synthetase_counts = synthetase_counts_by_group - (2 * synthetase_variance_by_group)
	scaled_synthetase_counts = synthetase_counts_by_group
	assert all(scaled_synthetase_counts > 0)

	predicted_trna_synthetase_rates = initialAAPolymerizationRate / scaled_synthetase_counts
	kb.trna_synthetase_rates = 2 * predicted_trna_synthetase_rates

	# fitKb2_metabolism(kb, simOutDir, bulkAverageContainer, bulkDeviationContainer)
