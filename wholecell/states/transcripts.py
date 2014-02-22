'''

Like chromosome.py, this is just a prototype.

'''

# Imports

from __future__ import division

import numpy
import tables


# Constants

MAX_TRANSCRIPTS = 2000 # expected max number of transcripts
TRANSCRIPT_LEN = 1000 # expected length of transcripts

ALLOCATION = 2 * MAX_TRANSCRIPTS * TRANSCRIPT_LEN # space allocated for transcripts

N_TRANSCRIPTS = 1000 # number of transcripts to create for the prototype simulation

RATE_TRANSCRIPTION = 42 # nt/s
N_SIM_STEPS = 10 # number of steps to simulate transcription


# Logging

eventLog = [] # Stores runtime events for later output

def printEventLog():
	# Collapses redundant events into single messages

	finalLog = []
	nCopies = 0
	for i, message in enumerate(eventLog[::-1]):
		index = len(eventLog) - i - 1

		if index > 0 and eventLog[index-1] == message:
			nCopies += 1

		else:
			if nCopies == 0:
				finalLog.append(message)

			else:
				finalLog.append(message + ' (x{})'.format(nCopies+1))

			nCopies = 0

	print ''
	print '\n'.join(finalLog[::-1])
	print ''


# Classes

class Transcript(object):
	# An mRNA molecule
	# TODO: reference this to a DNA sequence to determine the polymers' sequences and features (like promoters)
	extent = 0 # Amount transcribed

	startsAt = None # Starting location in Transcripts._array
	tscIndex = None # Location in Transcripts.transcripts list


class Transcripts(object):
	# TODO: track bound molecules, like the Chromosome object
	_empty = -2 # active regions with no bound molecules
	_inactive = -1 # inactive regions (allocated but not associated with any transcript)


	def __init__(self):
		self._array = numpy.empty(ALLOCATION, numpy.int64) # space allocated for transcripts to occupy
		self._array[:] = self._inactive

		self.transcripts = [] # existing transcripts

		eventLog.append('Initialized')


	def setup(self):
		# Create RNA transcripts of random length

		for i in xrange(N_TRANSCRIPTS):
			transcript = Transcript()
			transcript.extent = numpy.random.randint(TRANSCRIPT_LEN)
			self.transcriptAdd(transcript)

		eventLog.append('Transcripts set up')


	def transcriptAdd(self, transcript):
		# Add a new transcript

		tscIndex = self._transcriptAssignIndex(transcript)
		arrIndex = self._findFreeRegion(Transcript.extent)

		transcript.startsAt = arrIndex
		self._rangeIs(transcript.startsAt, transcript.startsAt+transcript.extent, self._empty)

		eventLog.append('Added a transcript')


	def transcriptExtend(self, transcript, length):
		# Extend the transcript, activating more regions for molecules to bind

		if not self._rangeInactive(
				transcript.startsAt + transcript.extent,
				transcript.startsAt + transcript.extent + length
				):
			
			eventLog.append('Relocating transcript due to collision on extension')

			self._rangeIs(transcript.startsAt, transcript.startsAt + transcript.extent, self._inactive)
			arrIndex = self._findFreeRegion(transcript.extent + length)
			transcript.startsAt = arrIndex


		transcript.extent += length
		self._rangeIs(transcript.startsAt, transcript.startsAt+transcript.extent, self._empty)


	def _transcriptAssignIndex(self, transcript):
		# Assign a new transcript to an index

		# NOTE: this logic is overly complex but should nicely handle removing
		# transcripts once that feature is desired.

		try:
			tscIndex = self.transcripts.index(None)

		except ValueError:
			self.transcripts.append(None)
			tscIndex = len(self.transcripts) - 1

		self.transcripts[tscIndex] = transcript
		transcript.tscIndex = tscIndex

		return tscIndex


	def _findFreeRegion(self, length):
		# Find a region of a desired length that is open

		# NOTE: I've intentionally made this sub-optimal to make it interesting. Even
		# so, after one round of elongation the number of collisions falls drastically,
		# suggesting that even "random" isn't a terrible heuristic.

		while True:
			arrIndex = numpy.random.randint(ALLOCATION)

			if self._rangeInactive(arrIndex, arrIndex + length):
				return arrIndex


	def _rangeInactive(self, start, stop):
		# Return whether a region on the array is fully inactive

		return (self._array[self._indexes(start, stop)] == self._inactive).all()


	def _rangeIs(self, start, stop, value):
		# Set a region on the array to a value

		self._array[self._indexes(start, stop)] = value


	def _indexes(self, start, stop):
		# Circularly-permuted indexing

		# TODO: determine overhead of this method

		if stop >= ALLOCATION:
			return numpy.r_[start:ALLOCATION, :stop%ALLOCATION]

		else:
			return numpy.s_[start:stop]


	def occupancy(self):
		# Calculate how full the allocation is (primarily for illustrative/debugging purposes)

		return 1 - (self._array == self._inactive).sum() / ALLOCATION


# Main execution

def main():
	# Build the state and run a few fake steps of transcriptions
	transcripts = Transcripts()
	transcripts.setup()

	# Extend every transcript
	for t in xrange(N_SIM_STEPS):
		eventLog.append('Simulating extension step #{}'.format(t))
		# TODO: initiation and termination of transcription, degradation of transcript

		for transcript in transcripts.transcripts: # I apologize for these horrible names
			transcripts.transcriptExtend(transcript, RATE_TRANSCRIPTION)

		# eventLog.append('Array occupancy is {:0.2%}'.format(transcripts.occupancy()))

	printEventLog()

	# import ipdb
	# ipdb.set_trace()


if __name__ == '__main__':
	main()
