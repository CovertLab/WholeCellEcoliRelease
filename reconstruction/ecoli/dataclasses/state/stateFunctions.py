import numpy as np
from wholecell.utils.unit_struct_array import UnitStructArray

def addToStateCommon(bulkState, ids, masses):
	newAddition = np.zeros(
		len(ids),
		dtype = [
			("id", "a50"),
			("mass", "{}f8".format(masses.asNumber().shape[1])), # TODO: Make this better
			]
		)

	bulkState.units['mass'].matchUnits(masses)

	newAddition["id"] = ids
	newAddition["mass"] = masses.asNumber()

	return UnitStructArray(np.hstack((bulkState.fullArray(), newAddition)), bulkState.units)

def createIdsInAllCompartments(ids, compartments):
	idsByCompartment = [
		'{}[{}]'.format(i, c)
		for c in compartments
		for i in ids
		]
	return np.array(idsByCompartment)

def createIdsWithCompartments(dictList):
	return ['{}[{}]'.format(x['id'], c)
			for x in dictList
			for c in x['location']
			]

def createMassesByCompartments(dictList):
	return np.array([x['mw']
			for x in dictList
			for c in x['location']
			], dtype = np.float64)