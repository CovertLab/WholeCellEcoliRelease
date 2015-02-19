import numpy as np

def addToBulkState(bulkState, ids, masses):
	newAddition = np.zeros(
		len(ids),
		dtype = [
			("id", "a50"),
			("mass", "{}f8".format(masses.shape[1])), # TODO: Make this better
			]
		)

	newAddition["id"] = ids
	newAddition["mass"] = masses
	return np.hstack((bulkState, newAddition))

def createIdsInAllCompartments(ids, compartments):
	idsByCompartment = [
		'{}[{}]'.format(i, c)
		for c in compartments
		for i in ids
		]
	return np.array(idsByCompartment)

