
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def starvationVariantTotalIndices(sim_data):
	nConditions = 8
	return nConditions

def starvationVariant(sim_data, index):
	# Vary code turning on and off for simulation

	if index == 0:
		return CONTROL_OUTPUT, sim_data
	elif index == 1:
		sim_data.translationSaturation = True
	elif index == 2:
		sim_data.glucoseLimitation = True
		sim_data.translationSaturation = True
		sim_data.fractionGlucoseLimit = 0.4
		sim_data.synthetase_km_scale = 0.8
	elif index == 3:
		sim_data.glucoseLimitation = True
		sim_data.translationSaturation = True
		sim_data.fractionGlucoseLimit = 0.6
		sim_data.synthetase_km_scale = 0.8
	elif index == 4:
		sim_data.glucoseLimitation = True
		sim_data.translationSaturation = True
		sim_data.fractionGlucoseLimit = 0.7
		sim_data.synthetase_km_scale = 0.8

	return dict(
		shortName = "glc={},transSat={}".format(sim_data.glucoseLimitation, sim_data.translationSaturation),
		desc = "glucoseLimitation={},by {},translationSaturation={}, f={}".format(sim_data.glucoseLimitation, sim_data.fractionGlucoseLimit, sim_data.translationSaturation, sim_data.synthetase_km_scale)
		), sim_data
