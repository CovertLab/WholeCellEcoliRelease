"""
This class contains dictionaries of common names and synonyms of simulation
elements in sim_data. Common names can be retrieved by calling
get_common_names() and a list of synonyms can be retrieved by calling
get_synonyms(). Both the common names and synonyms are gathered from relevant
columns in raw data files.
"""
import itertools

class CommonNames(object):
	def __init__(self, raw_data):
		self._build_common_names(raw_data)
		self._build_synonyms(raw_data)

	def get_common_name(self, element_id):
		"""
		Returns the common name of the simulation element with the given ID. If
		a common name does not exist, returns the original element ID.

		Args:
			element_id (str): ID of the simulation element

		Returns: Common name of the given simulation element (str)
		"""
		return self._common_names.get(element_id, element_id)

	def get_synonyms(self, element_id):
		"""
		Returns a list of synonyms for the simulation element with the given ID.
		If a list of synonyms does not exist, returns a singleton list with the
		original element ID as the only entry.

		Args:
			element_id (str): ID of the simulation element

		Returns: List of synonyms for the given simulation element (str)
		"""
		return self._synonyms.get(element_id, [element_id])


	def _build_common_names(self, raw_data):
		# Create dictionary of common names
		self._common_names = {}

		# Get common names of protein complexes
		for rxn in itertools.chain(
				raw_data.complexation_reactions,
				raw_data.equilibrium_reactions):
			if rxn['common_name'] is not None and len(rxn['common_name']) > 0:
				complex_name = None
				for (mol_name, stoich) in rxn['stoichiometry'].items():
					if stoich == 1:
						complex_name = mol_name

				if complex_name is not None:
					self._common_names[complex_name] = rxn['common_name']

		# Get common names (symbols) of genes
		for gene in raw_data.genes:
			if gene['symbol'] is not None and len(gene['symbol']) > 0:
				self._common_names[gene['id']] = gene['symbol']

		# Get common names of protein monomers, metabolites, RNAs, and transcription units
		for row in itertools.chain(
				raw_data.metabolites, raw_data.proteins, raw_data.rnas,
				raw_data.transcription_units):
			if row['common_name'] is not None and len(row['common_name']) > 0:
				self._common_names[row['id']] = row['common_name']

	def _build_synonyms(self, raw_data):
		# Create dictionary of synonyms
		self._synonyms = {}

		# Get synonyms of genes, metabolites, protein monomers, and RNAs
		for row in itertools.chain(
				raw_data.genes, raw_data.metabolites, raw_data.proteins, raw_data.rnas):
			if row['synonyms'] is not None and len(row['synonyms']) > 0:
				self._synonyms[row['id']] = row['synonyms']
