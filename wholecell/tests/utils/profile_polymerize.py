"""
profile_polymerize.py

Profiles where the time goes in polymerize().

Instructions: From the wcEcoli directory, run
	kernprof -lv wholecell/tests/utils/profile_polymerize.py
to get a line-by-line profile of functions decorated @profile.

This is split out into a separate file and directory from polymerize.py because
Python puts the current script's directory on the sys.path, and that breaks for
polymerize.py since its directory contains random.py. Importing numpy tries to
import the system random which instead loads the local random.py, and that
fails. So keep polymerize.py's directory off the sys.path.

The normal solution is to cd to wcEcoli then run "python -m <module>" but
kernprof doesn't support that.

@author: Jerry Morrison, John Mason, Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/10/2016
"""
from __future__ import absolute_import, division, print_function

import __builtin__
import sys
import os
import time
import cProfile
import pstats
import StringIO

import numpy as np

# EXPECTS: The current working directory is "wcEcoli/".
# Put that on the sys path in place of this script's directory so
# `import wholecell.utils.polymerize` and that module's imports will work even
# when run via kernprof.
sys.path[0] = os.getcwd()  # Put cwd on the python path so imports will work

from wholecell.utils.polymerize import polymerize

PAD_VALUE = polymerize.PAD_VALUE

# Wrap with kernprof profiling decorator - will throw an error if we call this
# script using the vanilla python interpreter.
if 'profile' not in __builtin__.__dict__:
	raise Exception(
		'kernprof @profile decorator not available.  This script should be '
		+ 'invoked via kernprof -lv.  If invoked correctly and this error '
		+ 'message is still raised, see issue #117.'
		)

# Attach __iter__ method to preserve old interface
# TODO (John): migrate to new interface
polymerize.__iter__ = lambda self: iter((self.sequenceElongation, self.monomerUsages, self.nReactions))

# Wrap methods in line-profiling decorator
# TODO (John): write an introspection utility to facilitate decoration

# noinspection PyUnresolvedReferences
def setup_profiler():
	polymerize.__init__ = profile(polymerize.__init__)

	polymerize._setup = profile(polymerize._setup)
	polymerize._sanitize_inputs = profile(polymerize._sanitize_inputs)
	polymerize._gather_input_dimensions = profile(polymerize._gather_input_dimensions)
	polymerize._gather_sequence_data = profile(polymerize._gather_sequence_data)
	polymerize._prepare_running_values = profile(polymerize._prepare_running_values)
	polymerize._prepare_outputs = profile(polymerize._prepare_outputs)

	polymerize._elongate = profile(polymerize._elongate)
	polymerize._elongate_to_limit = profile(polymerize._elongate_to_limit)
	polymerize._finalize_resource_limited_elongations = profile(polymerize._finalize_resource_limited_elongations)
	polymerize._update_elongation_resource_demands = profile(polymerize._update_elongation_resource_demands)

	polymerize._finalize = profile(polymerize._finalize)
	polymerize._clamp_elongation_to_sequence_length = profile(polymerize._clamp_elongation_to_sequence_length)

setup_profiler()

def _setupRealExample():
	# Test data pulled from an actual sim at an early time point.
	monomerLimits = np.array([
		11311,  6117,  4859,  6496,   843,  7460,  4431,  8986,  2126,
		6385,  9491,  7254,  2858,  3770,  4171,  5816,  6435,  1064,
		3127,     0,  8749])

	randomState = np.random.RandomState()

	nMonomers = len(monomerLimits) # number of distinct aa-tRNAs
	nSequences = 10000 # approximate number of ribosomes
	length = 16 # translation rate
	nTerminating = np.int64(length/300 * nSequences) # estimate for number of ribosomes terminating

	sequences = np.random.randint(nMonomers, size = (nSequences, length))

	sequenceLengths = length * np.ones(nSequences, np.int64)
	sequenceLengths[np.random.choice(nSequences, nTerminating, replace = False)] = np.random.randint(length, size = nTerminating)

	sequences[np.arange(length) > sequenceLengths[:, np.newaxis]] = PAD_VALUE

	reactionLimit = 10000000

	return sequences, monomerLimits, reactionLimit, randomState

def _setupExample():
	# Contrive a scenario which is similar to real conditions

	randomState = np.random.RandomState()

	nMonomers = 36 # number of distinct aa-tRNAs
	nSequences = 10000 # approximate number of ribosomes
	length = 16 # translation rate
	nTerminating = np.int64(length/300 * nSequences) # estimate for number of ribosomes terminating
	monomerSufficiency = 0.85
	energySufficiency = 0.85

	sequences = np.random.randint(nMonomers, size = (nSequences, length))

	sequenceLengths = length * np.ones(nSequences, np.int64)
	sequenceLengths[np.random.choice(nSequences, nTerminating, replace = False)] = np.random.randint(length, size = nTerminating)

	sequences[np.arange(length) > sequenceLengths[:, np.newaxis]] = PAD_VALUE

	maxReactions = sequenceLengths.sum()

	monomerLimits = monomerSufficiency * maxReactions/nMonomers*np.ones(nMonomers, np.int64)
	reactionLimit = energySufficiency * maxReactions

	return sequences, monomerLimits, reactionLimit, randomState


def _simpleProfile():
	np.random.seed(0)

	sequences, monomerLimits, reactionLimit, randomState = _setupExample()

	nSequences, length = sequences.shape
	nMonomers = monomerLimits.size
	sequenceLengths = (sequences != PAD_VALUE).sum(axis = 1)
	elongation_rates = []  # TODO: What to use here?

	t = time.time()
	sequenceElongation, monomerUsages, nReactions = polymerize(
		sequences, monomerLimits, reactionLimit, randomState, elongation_rates)
	evalTime = time.time() - t

	assert (sequenceElongation <= sequenceLengths+1).all()
	assert (monomerUsages <= monomerLimits).all()
	assert nReactions <= reactionLimit
	assert nReactions == monomerUsages.sum()

	print("""
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
		monomerUsages.sum()/monomerLimits.sum(),
		nReactions/reactionLimit,
		(sequenceElongation == sequenceLengths).sum()/nSequences,
		sequenceElongation.sum()/sequenceLengths.sum()
		))


def _fullProfile():
	np.random.seed(0)

	sequences, monomerLimits, reactionLimit, randomState = _setupRealExample()

	# Recipe from https://docs.python.org/2/library/profile.html#module-cProfile
	pr = cProfile.Profile()
	pr.enable()

	elongation_rates = []  # TODO: What to use here?
	sequenceElongation, monomerUsages, nReactions = polymerize(
		sequences, monomerLimits, reactionLimit, randomState, elongation_rates)

	pr.disable()
	s = StringIO.StringIO()
	sortby = 'cumulative'
	ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
	ps.print_stats()
	print(s.getvalue())


if __name__ == "__main__":
	_fullProfile()
