import copy
import numpy as np

trnaCapacity = np.random.randint(100, size = (21,))
synthetaseCapacity = np.random.randint(100, size = (21,))
aaCountInSequence = np.max(np.vstack((trnaCapacity, synthetaseCapacity)), axis = 0) + np.random.randint(50, size = (21,))
aaCountInSequence[0] = trnaCapacity[0] - 50
aaCounts = copy.copy(aaCountInSequence)


aaLimitation = (aaCountInSequence - aaCounts).clip(min = 0)
trnaCapacityLimitation = (aaCountInSequence - trnaCapacity).clip(min = 0)
synthetaseCapacityLimitation = (aaCountInSequence - synthetaseCapacity).clip(min = 0)

aaExcess = -1 * (aaCountInSequence - aaCounts).clip(max = 0)
trnaCapacityExcess = -1 * (aaCountInSequence - trnaCapacity).clip(max = 0)
synthetaseCapacityExcess = -1 * (aaCountInSequence - synthetaseCapacity).clip(max = 0)

# 0 - not stalled
# 1 - aa limitation
# 2 - trna limitation
# 3 - synthetase limitation

##############################################

trnaCapacity = np.random.randint(100, size = (50, 21))
synthetaseCapacity = np.random.randint(100, size = (50,21))
aaCountInSequence = np.max(np.vstack((trnaCapacity, synthetaseCapacity)), axis = 0) + np.random.randint(50, size = (50,21))
aaCountInSequence[:,0] = trnaCapacity[:,0] - 5
aaCounts = copy.copy(aaCountInSequence)

aaLimitation = (aaCountInSequence - aaCounts).clip(min = 0).sum(axis = 1)
trnaCapacityLimitation = (aaCountInSequence - trnaCapacity).clip(min = 0).sum(axis = 1)
synthetaseCapacityLimitation = (aaCountInSequence - synthetaseCapacity).clip(min = 0).sum(axis = 1)

aaExcess = -1 * (aaCountInSequence - aaCounts).clip(max = 0).sum(axis = 1)
trnaCapacityExcess = -1 * (aaCountInSequence - trnaCapacity).clip(max = 0).sum(axis = 1)
synthetaseCapacityExcess = -1 * (aaCountInSequence - synthetaseCapacity).clip(max = 0).sum(axis = 1)