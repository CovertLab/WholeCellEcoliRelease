#!/usr/bin/env python

import numpy as np

def PolymerizeMatrix(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase):
	'''
	polymerize
	Computes the maximum polymerization for the available amounts of bases
	for all active polymers using a greedy algorithm.

	1. Finds the most limited bases and the position at which they become limiting.
	2. Finds the energy limit of polymerization.
	3. Polymerizes all polymers up to, but not including, the limiting position.
	4. Polymerizes polymers at the limiting position according to base availability.
	5. Culls polymers that cannot make additional progress.
	6. Repeats 1-5 until no additional bases can be polymerized.

	@type sequences:	numpy matrix of str
	@param sequences:	[['a','g','c','t','c','g','c','a','t'],['t','g','c','g','c','g','g','a','t'],...etc.]
	@type baseAmounts:	ndarray of int
	@param baseAmounts:	Amounts of each base. i.e. [10,22,4,5]
	@type bases:	ndarray of str
	@param bases:	['1','2','3','4','5','6',...]
	@type basePadValue:	str
	@param basePadValue:	' '
	@type energy:	 int
	@param energy:	Amount of energy available
	@type energyCostPerBase:	int
	@param energyCostPerBase:	Amount of energy required to polymerize one base
	@rtype:	tuple of (ndarray of int, ndarray of int, ndarray of int, int, int)
	@return:	(progress = number of bases polymerized in each sequence, 
				baseAmounts = amounts of each base after polymerization,
				baseCosts = amount of each base used (should sum to input with baseAmounts),
				energy = amount of energy available after polymerization,
				energyCost = amount of energy used (should sum to input with energy))
	'''
	
	# Validate input
	validatePolymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase)
	
	nSequences		=	sequences.shape[0]						#number of sequences
	baseCosts		=	np.zeros((len(bases),), dtype = np.int)	#usage of bases, in order of bases
	energyCost		=	0										#usage of energy
	progress		=	np.zeros((nSequences,), dtype = np.int)	#number of bases polymerized in each sequence
	activeSeqIdxs	=	np.array(range(nSequences))				#indexs of sequences still being polymerized
	seqLengths		=	lengths(sequences, basePadValue)		#lengths of sequences
	
	
	# Elongate sequences
	while sequences.size != 0 and energy >= energyCostPerBase:
		# Eliminate sequences with nothing to polymerize
		toKeep_tf = progress[activeSeqIdxs] < seqLengths[activeSeqIdxs]
		activeSeqIdxs = activeSeqIdxs[toKeep_tf]
		sequences = sequences[toKeep_tf, :]
		if sequences.size == 0:
			# Break if no more sequences
			break
		
		# Calculate limit of elongation and the number of bases each of the active sequences can be elongated
		(elongation, baseUsage, limitingBases) = calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase)
	
		# Elongate all active sequences up to elongation
		# progress is the count of the number of bases that have been polymerized in each sequence
		progress[activeSeqIdxs] = np.min(np.vstack((progress[activeSeqIdxs] + elongation, seqLengths[activeSeqIdxs])), axis = 0)
		baseAmounts = baseAmounts - baseUsage
		baseCosts = baseCosts + baseUsage
		energy = energy - (energyCostPerBase * np.sum(baseUsage))
		energyCost = energyCost + (energyCostPerBase * np.sum(baseUsage))
		sequences = sequences[:, elongation:]
		
		if sequences.size == 0:
			# Break if no more sequences
			break
		
		if limitingBases.size != 0:
			# Eliminate sequences for which bases are not available
			toKeep_tf = np.array([True]*sequences.shape[0])
			for i in range(len(limitingBases)):
				# Indicies of reactions limited by base i
				seqIdxs = np.where(sequences[:,0] == bases[limitingBases[i]])[0]
				# Randomize order of indicies of reactions limited by base i
				randOrder = np.random.permutation(seqIdxs)
				# Keep first n reactions, where n is the number of the limiting bases remaining, and set the rest to False so
				# they are removed
				# Example:
				#	baseAmounts = 2
				#	randOrder = [2,4,1,3]
				#	Remove reactions with indicies 1 and 3 and keep 2 and 4. While loop will start over missing those
				#	two reactions and elongation could proceed farther.
				toKeep_tf[randOrder[baseAmounts[limitingBases[i]]:]] = False
			activeSeqIdxs = activeSeqIdxs[toKeep_tf]
			sequences = sequences[toKeep_tf, :]
		elif elongation == 0:
			# Limited energy and no limiting bases. Randomly select sequences to receive one last base
			n = int(np.fix(energy / energyCostPerBase))
			idxs = np.random.permutation(len(activeSeqIdxs))
			idxs_addTo = idxs[:n]
			progress[activeSeqIdxs[idxs_addTo]] = progress[activeSeqIdxs[idxs_addTo]] + 1
			baseUsage = countBases(sequences[idxs_addTo,0], bases)
			baseUsage = np.array(baseUsage[:,1].T)[0]
			baseAmounts = baseAmounts - baseUsage
			baseCosts = baseCosts + baseUsage
			energy = energy - energyCostPerBase * n
			energyCost = energyCost + energyCostPerBase * n
			break
		
	return progress, baseAmounts, baseCosts, energy, energyCost
	
	
def validatePolymerize(sequences, baseAmounts, bases, basePadValue, energy, energyCostPerBase):
	s = ''
	furtherValidate = True
	if not isinstance(sequences, np.matrix):
		s += 'sequences must be a numpy matrix!\n'
		furtherValidate = False
	elif sequences.dtype.type != np.string_ and sequences.dtype.type != np.unicode_:
		s += 'sequences must be matrix of strings!\n'
		furtherValidate = False
	
	if not isinstance(baseAmounts, np.ndarray):
		s += 'baseAmounts must be an ndarray!\n'
		furtherValidate = False
	elif baseAmounts.dtype != np.int:
		s += 'baseAmounts must be ndarray of ints!\n'
		furtherValidate = False
	elif len(baseAmounts.shape) != 1:
		s += 'baseAmounts must be one dimensional!\n'
		furtherValidate = False
	
	if not isinstance(bases, np.ndarray):
		s += 'bases must be an ndarray!\n'
		furtherValidate = False
	elif bases.dtype.type != np.string_:
		s += 'bases must be an ndarray of strings!\n'
		furtherValidate = False
	elif len(bases.shape) != 1:
		s += 'bases must be one dimensional!\n'
		furtherValidate = False
	
	if not isinstance(basePadValue, str):
		s += 'basePadValue must be a string!\n'
		furtherValidate = False
	
	if not isinstance(energy, int):
		s += 'energy must be an int!\n'
	elif energy < 0:
		s += 'energy must be greater than zero!\n'
	
	if not isinstance(energyCostPerBase, int):
		s += 'energyCostPerBase must be an int!\n'
	elif energyCostPerBase < 0:
		s += 'energyCostPerBase must be greater than zero!\n'
	
	if furtherValidate:
		# Validate seqeunces and alphabet
		if sequences.size == 0:
			s += 'sequences has no data!\n'
		
		usedAlphabet = [0]*sequences.shape[0]*sequences.shape[1]
		count = 0
		for row in range(sequences.shape[0]):
			for col in range(sequences.shape[1]):
				usedAlphabet[count] = sequences[row,col]
				count += 1
		usedAlphabet = "".join(usedAlphabet)
		lenUsed = len(usedAlphabet)
		charCount = 0
		for value in bases:
			charCount += usedAlphabet.count(value)
		charCount += usedAlphabet.count(basePadValue)
		
		if charCount < lenUsed:
			s += 'sequences uses incorrect alphabet!\n'
		
		# Validate basePadValue
		badPad = False
		if isinstance(bases, np.ndarray):
			for element in bases:
				if basePadValue == element:
					badPad = True
		else:
			for letter in bases:
				if basePadValue == letter:
					badPad = True
		if badPad:
			s += 'basePadValue has the same value as a base!\n'
	
	if len(s):
		raise polymerizeMatrixException, s


def calculateElongationLimits(sequences, bases, baseAmounts, energy, energyCostPerBase):
	'''
	calculateElongationLimits
	
	@type sequences:	numpy matrix of str
	@param sequences:	[['a','g','c','t','c','g','c','a','t'],['t','g','c','g','c','g','g','a','t'],...etc.]
	@type bases:	ndarray of str
	@param bases:	['1','2','3','4','5','6',...]
	@type baseAmounts:	ndarray
	@param baseAmounts:	Amounts of each base. i.e. [10,22,4,5]
	@type energy:	 int
	@param energy:	Amount of energy available
	@type energyCostPerBase:	int
	@param energyCostPerBase:	Amount of energy required to polymerize one base
	@rtype:	tuple of (int, ndarray, ndarray)
	@return:	(length of successful elongation, count of each base used, indicies of limting bases)
	'''
	
	# Calculate cumulative base counts across sequences
	cumBaseCounts = np.cumsum(countBases(sequences, bases), axis = 1)
	nBases = cumBaseCounts.shape[0]
	
	# (1) Individual base limits - first index where bases are limiting
	# polymerization could proceed for each base to base limit - 1 index
	baseElongationLimits = np.zeros((nBases,), dtype = np.int)
	for i in range(nBases):
		baseElongationLimits[i] = np.max(np.where(cumBaseCounts[i,:] <= baseAmounts[i])[1])
	
	# (2) Energy limit - first index where energy is limiting
	# polymerization could proceed to energy limit - 1 index
	
	# Builds list of columb indicis where energy is not limiting
	energyElongationLimit = np.max(np.where(np.sum(cumBaseCounts, axis = 0) * energyCostPerBase <= energy)[1])
	
	# (3) Elongate up to minimum of individual base and energy limits
	# elongation is last position that could be polymerized up until
	elongation = np.min(np.append(baseElongationLimits, np.array([energyElongationLimit])))
	baseUsage = np.array(cumBaseCounts[:, elongation].T)[0]
	
	# (4) If acceptable, find bases which are setting the elongation limit
	if elongation < sequences.shape[1]:
		limitingBases = np.where(baseElongationLimits == elongation)[0]
	else:
		limitingBases = np.array([])
	
	return elongation, baseUsage, limitingBases


def countBases(sequences, bases):
	'''
	countBases
	Counts the number of each base used at each position in all of the sequences
	
	@type sequences:	numpy matrix of str
	@param sequences:	[['a','g','c','t','c','g','c','a','t'],['t','g','c','g','c','g','g','a','t'],...etc.]
	@type bases:	ndarray of str
	@param bases:	['1','2','3','4','5','6',...]
	@rtype:	numpy matrix
	@return:	rows correspond to bases, columns to number of each base at that position in all sequences
	'''
	
	# NOTE: Count bases appends ane extra column of zeros on the LHS of the matrix. This
	# makes keeping track of elongation easier.
	counts = np.matrix(np.zeros((len(bases), sequences.shape[1] + 1), dtype = np.int))
	for i in range(len(bases)):
		counts[i,1:] = np.sum(sequences == bases[i], axis = 0)
	return counts

	
def lengths(sequences, padValue):
	'''
	lengths
	Calculates length of each sequence
	
	@type sequences:	numpy matrix of str
	@param sequences:	[['a','g','c','t','c','g','c','a','t'],['t','g','c','g','c','g','g','a','t'],...etc.]
	@type padValue:	str
	@param padValue:	' '
	@rtype:	ndarray
	@return:	lengths of each sequence
	'''
	
	present_tf = sequences != padValue
	n = sequences.shape[0]
	result = np.zeros((n,), dtype = np.int)
	for i in range(n):
		result[i] = np.max([0, len(np.where(present_tf[i,:])[0])])
	return result


class polymerizeMatrixException(Exception):
	"""
	polymerizeMatrixException
	"""

if __name__ == "__main__":
	import time

	np.random.seed(0)

	characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghij_"

	# Contrive a scenario which is similar to real conditions

	nMonomers = 36 # number of distinct aa-tRNAs
	nSequences = 10000 # approximate number of ribosomes
	length = 16 # translation rate
	nTerminating = np.int64(1.* length/300 * nSequences) # estimate for number of ribosomes terminating
	monomerSufficiency = 0.85
	energySufficiency = 0.85
	monomers = np.arange(nMonomers)
	padValue = -1
	costPerMonomer = 1

	sequences = np.random.randint(nMonomers, size = (nSequences, length))

	sequenceLengths = length * np.ones(nSequences, np.int64)
	sequenceLengths[np.random.choice(nSequences, nTerminating, replace = False)] = np.random.randint(length, size = nTerminating)

	sequences[np.arange(length) > sequenceLengths[:, np.newaxis]] = padValue

	maxReactions = sequenceLengths.sum()

	monomerLimits = (monomerSufficiency * maxReactions/nMonomers*np.ones(nMonomers)).astype(np.int64)
	reactionLimit = np.int64(energySufficiency * maxReactions)

	# Cast to types this function supports
	sequences = np.matrix([[characters[val] for val in row] for row in sequences])
	monomers = np.array([characters[val] for val in monomers])
	padValue = characters[padValue]

	t = time.time()
	sequenceElongation, _, monomerUsages, _, nReactions = PolymerizeMatrix(
		sequences,
		monomerLimits.copy(),
		monomers,
		padValue,
		reactionLimit,
		costPerMonomer
		)
	evalTime = time.time() - t

	assert (sequenceElongation <= sequenceLengths+1).all()
	assert (monomerUsages <= monomerLimits).all()
	assert nReactions <= reactionLimit
	assert nReactions == monomerUsages.sum()

	print """
Polymerize function report:

For {} sequences of {} different monomers elongating by at most {}:

{:0.1f} ms to evaluate
{} polymerization reactions
{:0.1f} average elongations per sequence
{:0.1%} monomer utilization
{:0.1%} energy utilization
{:0.1%} fully elongated
{:0.1%} completion
""".format(
		nSequences,
		nMonomers,
		length,
		evalTime * 1000,
		nReactions,
		sequenceElongation.mean(),
		1.*monomerUsages.sum()/monomerLimits.sum(),
		1.*nReactions/reactionLimit,
		1.*(sequenceElongation == sequenceLengths).sum()/nSequences,
		1.*sequenceElongation.sum()/sequenceLengths.sum()
		)
