
# import in simulation data for mass, molecule counts

# plot mass/fractions against time (relative or absolute?)

# plot free NTPs vs time?

import os

import tables
import numpy as np
import matplotlib.pyplot as plt

SIM_DIR = os.path.join('out', 'plotting')

NTP_IDS = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]']

SIM_LENGTH = 3600

START_AT = 0
STRIDE = 10

MASS_NAMES = ["cellDry", "rna", "rrna", "trna", "mrna"]

timePoints = np.arange(START_AT, SIM_LENGTH + 1, STRIDE)

nTimePoints = timePoints.size

ntpCounts = np.zeros((nTimePoints, len(NTP_IDS)), np.int)

assignedIndexes = False
ntpIndexes = None

# with tables.open_file(os.path.join(SIM_DIR, 'BulkMolecules.hdf')) as h5file:
# 	if not assignedIndexes:
# 		names = h5file.get_node('/names')

# 		moleculeIds = names.moleculeIDs.read()

# 		ntpIndexes = np.array([moleculeIds.index(ntpId) for ntpId in NTP_IDS], np.int)

# 		assignedIndexes = True

# 	bulkMolecules = h5file.root.BulkMolecules

# 	ntpCounts = bulkMolecules.read(START_AT, None, STRIDE, 'counts')[:, ntpIndexes]

# plt.figure()
# plt.plot(timePoints, ntpCounts)

masses = np.zeros((nTimePoints, len(MASS_NAMES)), np.float)

with tables.open_file(os.path.join(SIM_DIR, 'Mass.hdf')) as h5file:
	mass = h5file.root.Mass

	for i, massName in enumerate(MASS_NAMES):
		masses[:, i] = mass.read(START_AT, None, STRIDE, massName)

masses /= masses[0, :]

plt.figure()
plt.plot(timePoints, masses)
plt.legend(MASS_NAMES, 'best')

plt.savefig('mass_fractions.png')
