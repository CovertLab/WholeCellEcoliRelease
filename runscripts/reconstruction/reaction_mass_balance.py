
THRESHOLD = 1e-9

from reconstruction.ecoli.knowledge_base import KnowledgeBaseEcoli

import numpy as np

kb = KnowledgeBaseEcoli(deleteLoadingData = False)

S = kb.process.metabolism.stoichMatrix()

masses = kb.getter.getMass(kb.process.metabolism.moleculeNames)

reactionNetMass = np.dot(S.T, masses.asNumber())

reactionIds = kb.process.metabolism.reactionIds

reactionNames = kb.process.metabolism.reactionNames

sorting = np.argsort(np.abs(reactionNetMass))[::-1]

print "-"*79

moleculeMass = dict([(x['id'],x['mw7.2']) for x in kb._metabolites])
moleculeFormula = dict([(x['id'],x['formula7.2']) for x in kb._metabolites])

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

	stoich = [x['stoichiometry'] for x in kb._reactions if reactionIds[index] == x['id']][0]
	stoich_short = [(x['coeff'],x['molecule'],moleculeMass[x['molecule']],moleculeFormula[x['molecule']]) for x in stoich]

	print "{}\t{}\t{}\t{}".format(reactionIds[index], reactionNames[index], reactionNetMass[index], stoich_short)

	i += 1

print "-"*79
print "Summary:"
print "{} bad reactions (out of {})".format(i-1, len(reactionNames))
print "{} reactions were under the threshold of {}".format(nSkipped, THRESHOLD)
print "{} reactions contained the substring 'exchange' or 'sink'".format(nExchange)
print "-"*79
