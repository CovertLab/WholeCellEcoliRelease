from __future__ import absolute_import, division, print_function

import numpy as np
from wholecell.utils.unit_struct_array import UnitStructArray

def addToStateCommon(bulkState, ids, masses):
	masses_unitless = masses.asNumber()

	if masses_unitless.ndim == 1:
		assert len(ids) == 1
		mass_size = len(masses_unitless)
	else:
		mass_size = masses_unitless.shape[1]

	newAddition = np.zeros(
		len(ids),
		dtype = [
			("id", "U50"),
			("mass", "{}f8".format(mass_size)),
			]
		)

	bulkState.units['mass'].matchUnits(masses)

	newAddition["id"] = ids
	newAddition["mass"] = masses.asNumber()

	return UnitStructArray(np.hstack((bulkState.fullArray(), newAddition)), bulkState.units)
