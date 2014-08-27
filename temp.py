import copy
import numpy as np

trnaCapacity = np.random.randint(100, size = (21,))
synthetaseCapacity = np.random.randint(100, size = (21,))
aaCountInSequence = np.max(np.vstack((trnaCapacity, synthetaseCapacity)), axis = 0) + np.random.randint(50, size = (21,))
aaCounts = copy.copy(aaCountInSequence)


allLimits = np.vstack((aaCountInSequence, aaCounts, trnaCapacity, synthetaseCapacity))
limitIndex = np.argmin(allLimits, axis = 0)
limitValue = np.min(allLimits, axis = 0)

stalls = aaCountInSequence - limitValue

attributeStalls = np.zeros(4)
attributeStalls[limitIndex] = stalls
# 0 - not stalled
# 1 - aa limited
# 2 - trna limited
# 3 - synthetase limited

##############################################


trnaCapacity = np.random.randint(100, size = (50, 21))
synthetaseCapacity = np.random.randint(100, size = (50,21))
aaCountInSequence = np.max(np.vstack((trnaCapacity, synthetaseCapacity)), axis = 0) + np.random.randint(50, size = (50,21))
aaCounts = copy.copy(aaCountInSequence)

allLimits = np.dstack((aaCountInSequence, aaCounts, trnaCapacity, synthetaseCapacity))
limitIndex = np.argmin(allLimits, axis = 2)
limitValue = np.min(allLimits, axis = 2)

stalls = aaCountInSequence - limitValue
attributeStalls = np.zeros((50,4))
attributeStalls[limitIndex] = stalls