
"""

TODO:
- document math
- replace fake metabolite pools with measured metabolite pools
- raise/warn if physiological metabolite pools appear to be smaller than what
 is needed at this time step size

"""

from __future__ import division

from itertools import izip

import numpy as np

from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from reconstruction.ecoli.fitter import countsFromMassAndExpression
from reconstruction.ecoli.fitter import normalize

def calcInitialConditions(sim, kb):
	randomState = sim.randomState

	timeStep = sim.timeStepSec() # This is a poor solution but will suffice for now

	bulkMolCntr = sim.states['BulkMolecules'].container

	# Set up states
	initializeBulkMolecules(bulkMolCntr, kb, randomState, timeStep)


def initializeBulkMolecules(bulkMolCntr, kb, randomState, timeStep):

	## Set other biomass components
	initializeBulkComponents(bulkMolCntr, kb, randomState, timeStep)


def initializeBulkComponents(bulkMolCntr, kb, randomState, timeStep):

	massFractions = kb.cellDryMassComposition[
		kb.cellDryMassComposition["doublingTime"].to("minute").magnitude == 60.0
		].fullArray()

	initDryMass = kb.avgCellDryMassInit.to("DCW_gram").magnitude
	cellMass = (
		kb.avgCellDryMassInit.to("DCW_gram").magnitude
		# + kb.avgCellWaterMassInit.magnitude
		)

	poolIds = kb.metabolitePoolIDs[:]

	mass = initDryMass
	mass -= massFractions["glycogenMassFraction"] * initDryMass
	mass -= massFractions["mureinMassFraction"] * initDryMass
	mass -= massFractions["lpsMassFraction"] * initDryMass
	mass -= massFractions["lipidMassFraction"] * initDryMass
	mass -= massFractions["inorganicIonMassFraction"] * initDryMass
	mass -= massFractions["solublePoolMassFraction"] * initDryMass

	# We have to remove things with zero concentration because taking the inverse of zero isn't so nice.
	poolIds = [x for idx, x in enumerate(kb.metabolitePoolIDs) if kb.metabolitePoolConcentrations.magnitude[idx] > 0]
	poolConcentrations = np.array([x for x in kb.metabolitePoolConcentrations.magnitude if x > 0])

	cellVolume = cellMass / kb.cellDensity
	cellDensity = kb.cellDensity.to("g / L").magnitude
	mws = kb.getMass(poolIds).to("g / mol").magnitude
	concentrations = poolConcentrations.copy()

	diag = cellDensity / (mws * concentrations) - 1
	A = -1 * np.ones((diag.size, diag.size))
	A[np.diag_indices(diag.size)] = diag
	b = mass * np.ones(diag.size)


	massesToAdd = np.linalg.solve(A, b)
	countsToAdd = massesToAdd / mws * kb.nAvogadro.to("1 / mol").magnitude

	V = (mass + massesToAdd.sum()) / cellDensity

	assert np.allclose(countsToAdd / kb.nAvogadro.magnitude / V, poolConcentrations)

	bulkMolCntr.countsIs(
		countsToAdd,
		poolIds
		)

	subunits = bulkMolCntr.countsView(["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"])
	subunitStoich = np.array([1, 1, 1])
	activeRibosomeMax = (subunits.counts() // subunitStoich).min()
	elngRate = kb.ribosomeElongationRate.to('amino_acid / s').magnitude
	T_d = kb.cellCycleLen.to("s").magnitude
	dt = kb.timeStep.to("s").magnitude

	activeRibosomesLastTimeStep = activeRibosomeMax * np.exp( np.log(2) / T_d * (T_d - dt)) / 2
	gtpsHydrolyzedLastTimeStep = activeRibosomesLastTimeStep * elngRate * kb.gtpPerTranslation

	bulkMolCntr.countsInc(
		gtpsHydrolyzedLastTimeStep,
		["GDP[c]"]
		)
