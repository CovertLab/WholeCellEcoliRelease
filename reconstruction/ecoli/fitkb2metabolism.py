import tables
from wholecell.utils.modular_fba import FluxBalanceAnalysis
def fitKb_2_metabolism(kb, simOutDir, bulkAverageContainer, bulkDeviationContainer):

	# Load the simulation output

	## Effective biomass reaction
	with tables.open_file(os.path.join(simOutDir, "ConcentrationChange.hdf")) as h5file:
		time = h5file.root.ConcentrationChange.col("time")
		timeStep = h5file.root.ConcentrationChange.col("timeStep")

		# NOTE: units are M/s
		concentrationChange = h5file.root.ConcentrationChange.col("concentrationChange")

		names = h5file.root.names
		biomassMoleculeIDs = np.array(names.moleculeIDs.read())

	## Find the most extreme concentration flux, after removing the first few time steps

	# TODO: intelligent filtering - use the correlation coefficient?

	concentrationChange = concentrationChange[3:, :] # NOTE: magic number

	concentrationChangeMostPositive = concentrationChange.max(0)
	concentrationChangeMostNegative = concentrationChange.min(0)

	effectiveBiomassReaction = concentrationChangeMostPositive.copy()

	negativeIsMostExtreme = (np.abs(concentrationChangeMostNegative)
		> concentrationChangeMostPositive)

	effectiveBiomassReaction[negativeIsMostExtreme] = concentrationChangeMostNegative[negativeIsMostExtreme]

	## Build the enzyme-fitting problem

	# TODO: write a class for setting up LP problems

	values = []
	rowIndexes = []
	colIndexes = []

	lowerValues = []
	lowerIndexes = []

	upperValues = []
	upperIndexes = []

	objValues = []
	objIndexes = []

	rowNames = []
	colNames = []

	### Set up reverse reactions

	dt = 1 * units.s

	reactionStoich = kb.process.metabolism.reactionStoich.copy()
	reactionEnzymes = kb.process.metabolism.reactionEnzymes.copy()
	reactionRates = kb.process.metabolism.reactionRates(dt)

	for reactionID in kb.process.metabolism.reversibleReactions:
		reverseReactionID = "{} (reverse)".format(reactionID)
		assert reverseReactionID not in reactionStoich.viewkeys()
		reactionStoich[reverseReactionID] = {
			moleculeID:-coeff
			for moleculeID, coeff in reactionStoich[reactionID].viewitems()
			}

		if reactionID in reactionEnzymes.viewkeys():
			reactionEnzymes[reverseReactionID] = reactionEnzymes[reactionID]

		if reactionID in reactionRates.viewkeys():
			reactionRates[reverseReactionID] = reactionRates[reactionID]

	### Set up metabolites and biochemical reactions

	for reactionID, stoich in reactionStoich.viewitems():
		assert reactionID not in colNames
		reactionIndex = len(colNames)
		colNames.append(reactionID)

		for moleculeID, coeff in stoich.viewitems():
			try:
				moleculeIndex = rowNames.index(moleculeID)

			except ValueError:
				moleculeIndex = len(rowNames)
				rowNames.append(moleculeID)

			rowIndexes.append(moleculeIndex)
			colIndexes.append(reactionIndex)
			values.append(coeff)

	### Set up exchange reactions

	initWaterMass = kb.constants.avgCellWaterMassInit
	initDryMass = kb.constants.avgCellDryMassInit

	initCellMass = initWaterMass + initDryMass

	initCellVolume = initCellMass / kb.constants.cellDensity

	coefficient = initDryMass / initCellVolume * dt

	exchangeConstraints = kb.process.metabolism.exchangeConstraints(
		kb.process.metabolism.externalExchangeMolecules,
		coefficient,
		units.mmol / units.L
		)

	for moleculeID, constraint in zip(kb.process.metabolism.externalExchangeMolecules, exchangeConstraints):
		exchangeID = "{} exchange".format(moleculeID)

		assert exchangeID not in colNames
		exchangeIndex = len(colNames)
		colNames.append(exchangeID)

		moleculeIndex = rowNames.index(moleculeID)

		rowIndexes.append(moleculeIndex)
		colIndexes.append(exchangeIndex)
		values.append(-1.0)

		lowerIndexes.append(exchangeIndex)
		lowerValues.append(-min(constraint, 1e6))

	### Set up biomass reaction

	biomassID = "biomass reaction"
	assert biomassID not in colNames
	biomassIndex = len(colNames)
	colNames.append(biomassID)

	effectiveBiomassReaction *= 10**3 # change to mmol

	for moleculeID, coeff in zip(biomassMoleculeIDs, effectiveBiomassReaction):
		moleculeIndex = rowNames.index(moleculeID)

		rowIndexes.append(moleculeIndex)
		colIndexes.append(biomassIndex)
		values.append(-coeff)

		lowerIndexes.append(biomassIndex)
		lowerValues.append(+1) # must be capable of producing 100% of the biomass in a step

	### Set up enzyme usage

	enzymeRatesAll = collections.defaultdict(set)

	for reactionID, enzymeID in reactionEnzymes.viewitems():
		reactionRate = reactionRates[reactionID]

		enzymeRatesAll[enzymeID].add(reactionRate)

	enzymeIDs = enzymeRatesAll.keys()
	perEnzymeRates = {
		enzymeID:max(enzymeRates)
		for enzymeID, enzymeRates in enzymeRatesAll.viewitems()
		}

	minimalEnzymeCounts = np.fmax(
		bulkAverageContainer.counts(enzymeIDs) - 2 * bulkDeviationContainer.counts(enzymeIDs),
		0
		)

	enzymeConc = (
		1 / kb.constants.nAvogadro / initCellVolume * minimalEnzymeCounts
		).asNumber(units.mmol / units.L)

	fullEnzymeRates = {
		enzymeID:perEnzymeRates[enzymeID] * enzymeConc[index]
		for index, enzymeID in enumerate(enzymeIDs)
		}

	for enzymeID, rateConstraint in fullEnzymeRates.viewitems():
		assert enzymeID not in rowNames
		enzymeIndex = len(rowNames)
		rowNames.append(enzymeID)

		constraintID = "{} constraint".format(enzymeID)
		assert constraintID not in colNames
		constraintIndex = len(colNames)
		colNames.append(constraintID)

		excessID = "{} excess capacity".format(enzymeID)
		assert excessID not in colNames
		excessIndex = len(colNames)
		colNames.append(excessID)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(constraintIndex)
		values.append(+1.0)

		upperIndexes.append(constraintIndex)
		upperValues.append(rateConstraint)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(excessIndex)
		values.append(+1.0)

		objIndexes.append(excessIndex)
		objValues.append(+1.0) # TODO: weighting

	for reactionID, enzymeID in reactionEnzymes.viewitems():
		if reactionID not in reactionRates.viewkeys():
			raise Exception("This code was not designed to handle enzymatic constraints without annotated rates.")

		reactionIndex = colNames.index(reactionID)
		enzymeIndex = rowNames.index(enzymeID)

		rowIndexes.append(enzymeIndex)
		colIndexes.append(reactionIndex)
		values.append(-1)

	import cvxopt
	import cvxopt.solvers

	nRows = max(rowIndexes) + 1
	nCols = max(colIndexes) + 1

	assert len(values) == len(rowIndexes) == len(colIndexes)

	A = cvxopt.spmatrix(values, rowIndexes, colIndexes)

	b = cvxopt.matrix(np.zeros(nRows, np.float64))

	assert len(objIndexes) == len(objValues)

	objectiveFunction = np.zeros(nCols, np.float64)
	objectiveFunction[objIndexes] = objValues

	f = cvxopt.matrix(objectiveFunction)

	G = cvxopt.matrix(np.concatenate(
		[np.identity(nCols, np.float64), -np.identity(nCols, np.float64)]
		))

	assert len(upperIndexes) == len(upperValues)

	upperBound = np.empty(nCols, np.float64)
	upperBound.fill(1e6)
	upperBound[upperIndexes] = upperValues

	assert len(lowerIndexes) == len(lowerValues)

	lowerBound = np.empty(nCols, np.float64)
	lowerBound.fill(0)
	lowerBound[lowerIndexes] = lowerValues

	h = cvxopt.matrix(np.concatenate([upperBound, -lowerBound], axis = 0))

	oldOptions = cvxopt.solvers.options.copy()

	cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

	solution = cvxopt.solvers.lp(f, G, h, A, b, solver = "glpk")

	cvxopt.solvers.options.update(oldOptions)
