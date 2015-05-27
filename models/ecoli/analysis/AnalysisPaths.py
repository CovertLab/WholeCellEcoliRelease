from os import listdir
from os.path import isdir, join

class AnalysisPaths(object):
	def __init__(self, seed_dir):
		self._generation_data = self._get_generations(seed_dir)
		#self._lineage_tree = self._get_lineage(seed_dir)

	def _get_generations(self, directory):
		generation_files = [join(directory,f) for f in listdir(directory) if isdir(join(directory,f)) and "generation" in f]
		self.n_generations = len(generation_files)
		generations = [None] * len(generation_files)
		for gen_file in generation_files:
			generations[int(gen_file[gen_file.rfind('_') + 1:])] = self._get_individuals(gen_file)
		return generations

	def _get_individuals(self, directory):
		individual_files = [join(directory,f) for f in listdir(directory) if isdir(join(directory,f))]
		individuals = [None] * len(individual_files)
		for ind_file in individual_files:
			individuals[int(ind_file[ind_file.rfind('/')+1:])] = ind_file
		return individuals

	def _get_lineage(self, directory):
		lineage = self._make_node_recursive(0, 0, None)
		return lineage

	def _make_node_recursive(self, gen, n, parent):
		if gen > len(self._generation_data) - 1:
			return None
		else:
			newNode = node(
				generation = gen,
				directory = self._generation_data[gen][n],
				)
			newNode.parent = parent
			newNode.daughter1 = self._make_node_recursive(gen+1, 2*n, newNode)
			newNode.daughter2 = self._make_node_recursive(gen+1, 2*n+1, newNode)
			return newNode

	def getGeneration(self, gen):
		return self._generation_data[gen]

	def getAll(self):
		allDir = []
		for x in self._generation_data:
			allDir.extend(x)
		return allDir

class node(object):
	def __init__(self, generation, directory, parent  = None, daughter1 = None, daughter2 = None):
		self.generation = generation
		self.directory = directory
		self.parent = parent
		self.daughter1 = daughter1
		self.daughter2 = daughter2
