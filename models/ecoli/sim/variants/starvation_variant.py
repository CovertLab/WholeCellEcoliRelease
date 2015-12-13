
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def starvationVariantTotalIndices(kb):
	nConditions = 8
	return nConditions

def starvationVariant(kb, index):
	# Vary code turning on and off for simulation

	if index == 0:
		return CONTROL_OUTPUT
	elif index == 1:
		kb.translationSaturation = True
	elif index == 2:
		kb.glucoseLimitation = True
		kb.translationSaturation = True
		kb.fractionGlucoseLimit = 0.35
		kb.synthetase_km_scale = 0.7
	elif index == 3:
		kb.glucoseLimitation = True
		kb.translationSaturation = True
		kb.fractionGlucoseLimit = 0.4
		kb.synthetase_km_scale = 0.7
	elif index == 4:
		kb.glucoseLimitation = True
		kb.translationSaturation = True
		kb.fractionGlucoseLimit = 0.45
		kb.synthetase_km_scale = 0.7
	elif index == 5:
		kb.glucoseLimitation = True
		kb.translationSaturation = True
		kb.fractionGlucoseLimit = 0.35
		kb.synthetase_km_scale = 0.8
	elif index == 6:
		kb.glucoseLimitation = True
		kb.translationSaturation = True
		kb.fractionGlucoseLimit = 0.4
		kb.synthetase_km_scale = 0.8
	elif index == 7:
		kb.glucoseLimitation = True
		kb.translationSaturation = True
		kb.fractionGlucoseLimit = 0.45
		kb.synthetase_km_scale = 0.8

	return dict(
		shortName = "glc={},transSat={}".format(kb.glucoseLimitation, kb.translationSaturation),
		desc = "glucoseLimitation={},by {},translationSaturation={}, f={}".format(kb.glucoseLimitation, kb.fractionGlucoseLimit, kb.translationSaturation, kb.synthetase_km_scale)
		)
