
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


	def attrs(self, moleculeName, attributes, timepoints = None):
		if timepoints is None:
			timepoints = self._timepoints

		# Load from the h5file

		moleculeTable = getattr(self._h5file.root, moleculeName)

		# Find the active entries

		isActive = (moleculeTable.col("_entryState") == 1)
		times = moleculeTable.col("_time")[isActive]

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


if __name__ == '__main__':
	filePath = "out/working/UniqueMolecules.hdf"

	with UniqueMoleculesData(filePath) as data:
		timepoints = data.timepoints()

		counts = data.counts("activeRnaPoly")

		rnaIndexes = data.attr("activeRnaPoly", "rnaIndex")

		required, assigned = data.attrs("activeRnaPoly", ["requiredACGU", "assignedACGU"])

	import matplotlib.pyplot as plt

	plt.plot(timepoints, counts)
	plt.show()
