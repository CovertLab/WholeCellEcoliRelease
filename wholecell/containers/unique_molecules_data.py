
from __future__ import division

from itertools import izip

import tables
import numpy as np

# TODO: caching for the timepoint:entries mappings
# TODO: tests
# TODO: refactor this class definition to focus on the container, not the state

# NOTE: I want to keep this class as small as possible - any common operations
# like averaging, summing, etc. can't be optimized any better here than they
# could be in the analysis script, and are more likely to break if placed in
# this file

class UniqueMoleculesData(object):
	def __init__(self, filePath):
		self._h5file = tables.open_file(filePath, 'r')

		# Load the table names

		globalReferenceTable = self._h5file.root._globalReference

		self._moleculeNames = globalReferenceTable.attrs.collectionNames

		# Record valid timepoints
		times = globalReferenceTable.col("_time")

		self._timepoints = np.unique(times)


	def timepoints(self):
		# Return the valid timepoints
		return self._timepoints.copy()


	def counts(self, moleculeName):
		# Return the number of active entries

		moleculeTable = getattr(self._h5file.root, moleculeName)

		isActive = (moleculeTable.col("_entryState") == 1)
		times = moleculeTable.col("_time")[isActive]

		return np.array([
			(times == t).sum() for t in self._timepoints
			])


	def attrsByMolecule(self, moleculeName, attributes):
		# Return a list of numpy struct arrays, where each array is the history
		# of a specific molecule

		# TODO: exclude inactive?

		# Load from the h5file

		moleculeTable = getattr(self._h5file.root, moleculeName)
		activeIndexes = np.where(moleculeTable.col("_entryState") == 1)[0]

		# Collect the unique molecule IDs
		uniqueIdsCol = moleculeTable.read_coordinates(activeIndexes, "_uniqueId")
		uniqueIds = np.unique(uniqueIdsCol)

		# allFields = moleculeTable.read_coordinates(activeIndexes)[attributes]

		n = uniqueIds.size

		# Yield the output
		for i, uniqueId in enumerate(uniqueIds):
			# if i % 1000 == 0:
			# 	print 100*i/n

			entriesForMolecule = (uniqueIdsCol == uniqueId)

			entries = moleculeTable.read_coordinates(activeIndexes[entriesForMolecule])

			yield entries[attributes]

			# entries = allFields[entriesForMolecule]

			# yield entries


	def attrs(self, moleculeName, attributes, timepoints = None):
		if timepoints is None:
			timepoints = self._timepoints

		# Load from the h5file

		moleculeTable = getattr(self._h5file.root, moleculeName)

		# Find the active entries

		isActive = (moleculeTable.col("_entryState") == 1)
		times = moleculeTable.col("_time")[isActive]

		# TODO: invert the for-loop order so this doesn't need to be cached
		entriesForTimepoints = [
			(times == time) for time in timepoints
			]

		# Collect the values for each attribute
		output = []

		for attribute in attributes:
			entries = moleculeTable.col(attribute)[isActive]

			values = []
			
			# Collect the attribute values for each time point
			for entriesForTimepoint in entriesForTimepoints:
				values.append(entries[entriesForTimepoint])

			output.append(values)
		
		return output


	def attr(self, moleculeName, attribute, timepoints = None):
		return self.attrs(moleculeName, (attribute,), timepoints)[0]


	def extractTimepoint(self, timepoint):
		# Return a UniqueObjectsContainer for a specific time point
		raise NotImplementedError("Waiting to implement this until some refactoring")


	def close(self):
		# Explictly close the file, just in case that is something you want to do
		self._h5file.close()


	# Methods for use as a context manager
	def __enter__(self):
		return self

	def __exit__(self, *excinfo):
		self.close()


def bundleByFieldValue(structArray, fieldName, outputFields = None):
	assert outputFields is None or isinstance(outputFields, list) # strict type requirements due to numpy
	# Collect the unique values for the field

	fieldValues = structArray[fieldName]

	uniqueFieldValues = np.unique(structArray[fieldName])

	# Yield the output
	for uniqueFieldValue in uniqueFieldValues:
		entriesForMolecule = (fieldValues == uniqueFieldValue)

		entries = structArray[entriesForMolecule][outputFields]

		yield uniqueFieldValue, entries
