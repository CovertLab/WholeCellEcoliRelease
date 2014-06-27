
THRESHOLD = 1e-9

from wholecell.reconstruction.knowledge_base_ecoli import KnowledgeBaseEcoli

import numpy as np

kb = KnowledgeBaseEcoli()

S = kb.metabolismStoichMatrix()

masses = kb.getMass(kb.metabolismMoleculeNames)

reactionNetMass = np.dot(S.T, masses)

reactionIds = kb.metabolismReactionIds

reactionNames = kb.metabolismReactionNames

sorting = np.argsort(np.abs(reactionNetMass))[::-1]

print "-"*79

i = 0
nSkipped = 0
nExchange = 0
for index in sorting:
	if np.abs(reactionNetMass[index]) < THRESHOLD:
		nSkipped += 1
		continue

	if "exchange" in reactionNames[index] or "sink" in reactionNames[index].lower():
		nExchange += 1
		continue

	print "{}\t{}\t{}".format(reactionIds[index], reactionNames[index], reactionNetMass[index])

	i += 1

print "-"*79
print "Summary:"
print "{} bad reactions (out of {})".format(i-1, len(reactionNames))
print "{} reactions were under the threshold of {}".format(nSkipped, THRESHOLD)
print "{} reactions contained the substring 'exchange' or 'sink'".format(nExchange)
print "-"*79
