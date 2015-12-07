
CONTROL_OUTPUT = dict(
	shortName = "control",
	desc = "Control simulation"
	)

def starvationVariantTotalIndices(kb):
	nConditions = 4
	return nConditions

def starvationVariant(kb, index):
	# Vary code turning on and off for simulation

	if index == 0:
		return CONTROL_OUTPUT
	elif index == 1:
		kb.glucoseLimitation = True
	elif index == 2:
		kb.translationSaturation = True
	elif index == 3:
		kb.glucoseLimitation = True
		kb.translationSaturation = True

	return dict(
		shortName = "glc={},transSat={}".format(kb.glucoseLimitation, kb.translationSaturation),
		desc = "glucoseLimitation={},translationSaturation={}".format(kb.glucoseLimitation, kb.translationSaturation)
		)
