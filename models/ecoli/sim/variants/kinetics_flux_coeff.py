
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def kineticsFluxCoeffTotalIndices(sim_data):
	nGenes = sim_data.process.transcription.rnaData.fullArray().size
	nConditions = nGenes + 1
	return nConditions


def kineticsFluxCoeff(sim_data, index):
	# Sets the coefficient of the upper limit of flux used for enzymeKinetics estimates

	nConditions = kineticsFluxCoeffTotalIndices(sim_data)

	if index % nConditions == 0:
		return CONTROL_OUTPUT, sim_data

	coeffIndex = (index - 1) % nConditions

	fluxLimit = 100000000. * 2**(-index)

	sim_data.constants.kineticRateLimitFactorUpper = fluxLimit

	return dict(
		shortName = "upper_flux %.1f" % (fluxLimit),
		desc = "Rate limit upper coefficient set to %.2f." % (fluxLimit)
		), sim_data