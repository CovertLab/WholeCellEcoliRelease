"""
SimulationData moleculeGroups

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

from reconstruction.ecoli.dataclasses.state.stateFunctions import createIdsWithCompartments

class MoleculeGroups(object):
	"""
	Helper class to extract molecule IDs of "special" groups of molecules. All
	values returned are lists of strings.
	"""

	def __init__(self, raw_data, sim_data):
		self._buildMoleculeGroups(raw_data, sim_data)

	def _buildMoleculeGroups(self, raw_data, sim_data):
		moleculeGroups = {
			'ntpIds': ["ATP[c]","CTP[c]","GTP[c]","UTP[c]"],
			'dNtpIds': ["DATP[c]", "DCTP[c]", "DGTP[c]", "TTP[c]"],
			'rnapIds': ['RPOB-MONOMER[c]', 'RPOC-MONOMER[c]', 'EG10893-MONOMER[c]'],

			'polymerizedAA_IDs': [x['id'] for x in raw_data.polymerized
				if x['is_aa'] and not x['is_end']], # TODO: end weight
			'polymerizedNT_IDs': [x['id'] for x in raw_data.polymerized
				if x['is_ntp'] and not x['is_end']], # TODO: end weight
			'polymerizedDNT_IDs': [x['id'] for x in raw_data.polymerized
				if x['is_dntp'] and not x['is_end']], # TODO: end weight

			's30_proteins':	['EG10912-MONOMER[c]', 'EG10916-MONOMER[c]',
				'EG10906-MONOMER[c]', 'EG10914-MONOMER[c]', 'EG10909-MONOMER[c]',
				'EG10903-MONOMER[c]', 'EG10911-MONOMER[c]', 'EG10904-MONOMER[c]',
				'EG10900-MONOMER[c]', 'EG10901-MONOMER[c]', 'EG10905-MONOMER[c]',
				'EG10915-MONOMER[c]', 'EG10918-MONOMER[c]', 'EG10919-MONOMER[c]',
				'EG10907-MONOMER[c]', 'EG11508-MONOMER[c]', 'EG10908-MONOMER[c]',
				'EG10920-MONOMER[c]', 'EG10910-MONOMER[c]', 'EG10902-MONOMER[c]',
				'EG10917-MONOMER[c]', 'EG10913-MONOMER[c]'],
			's30_16sRRNA': ['RRSA-RRNA[c]', 'RRSB-RRNA[c]', 'RRSC-RRNA[c]',
				'RRSD-RRNA[c]', 'RRSE-RRNA[c]', 'RRSG-RRNA[c]', 'RRSH-RRNA[c]'],

			's50_proteinComplexes': ['CPLX0-3956[c]'],
			's50_proteins':	['EG10872-MONOMER[c]', 'EG10879-MONOMER[c]',
				'EG11232-MONOMER[c]', 'EG10877-MONOMER[c]', 'EG10876-MONOMER[c]',
				'EG10892-MONOMER[c]', 'EG10874-MONOMER[c]', 'EG50001-MONOMER[c]',
				'EG10875-MONOMER[c]', 'EG10884-MONOMER[c]', 'EG11231-MONOMER[c]',
				'EG10887-MONOMER[c]', 'EG10878-MONOMER[c]', 'EG10886-MONOMER[c]',
				'EG10870-MONOMER[c]', 'EG10889-MONOMER[c]', 'EG10891-MONOMER[c]',
				'EG10888-MONOMER[c]', 'EG50002-MONOMER[c]', 'EG10869-MONOMER[c]',
				'EG10882-MONOMER[c]', 'EG10883-MONOMER[c]', 'EG10885-MONOMER[c]',
				'EG10890-MONOMER[c]', 'EG10864-MONOMER[c]', 'EG10881-MONOMER[c]',
				'EG10865-MONOMER[c]', 'EG10868-MONOMER[c]', 'EG10880-MONOMER[c]',
				'EG10867-MONOMER[c]', 'EG10866-MONOMER[c]'],
			's50_23sRRNA': ['RRLA-RRNA[c]', 'RRLB-RRNA[c]', 'RRLC-RRNA[c]',
				'RRLD-RRNA[c]', 'RRLE-RRNA[c]', 'RRLG-RRNA[c]', 'RRLH-RRNA[c]'],
			's50_5sRRNA': ['RRFA-RRNA[c]', 'RRFB-RRNA[c]', 'RRFC-RRNA[c]',
				'RRFD-RRNA[c]', 'RRFE-RRNA[c]', 'RRFG-RRNA[c]', 'RRFH-RRNA[c]'],

			# TODO: 'EG10245-MONOMER[c]' (DNAP III subunit tau) should be added
			# to the list of trimer subunits once frame-shifting proteins are
			# produced.
			'replisome_trimer_subunits': ['CPLX0-2361[c]', 'CPLX0-3761[c]'],
			'replisome_monomer_subunits': ['CPLX0-3621[c]', 'EG10239-MONOMER[c]',
				'EG11500-MONOMER[c]', 'EG11412-MONOMER[c]'],

			'aaIDs': sim_data.amino_acid_1_to_3_ordered.values(),
			'fragmentNT_IDs': [x['id'].replace('Polymerized','Fragment')
				for x in raw_data.polymerized if x['is_ntp'] and not x['is_end']],

			'exoRnaseIds': ['EG11620-MONOMER[c]', 'G7175-MONOMER[c]',
				'EG10858-MONOMER[c]', 'EG10863-MONOMER[c]', 'EG11259-MONOMER[c]',
				'EG11547-MONOMER[c]', 'EG10746-MONOMER[c]', 'G7842-MONOMER[c]',
				'EG10743-MONOMER[c]'],
			'endoRnase_RnaIDs':	['EG10856_RNA[c]', 'EG10857_RNA[c]',
				'EG10859_RNA[c]', 'EG10860_RNA[c]', 'EG10861_RNA[c]',
				'EG10862_RNA[c]', 'EG11299_RNA[c]', 'G7175_RNA[c]',
				'G7365_RNA[c]'],
			'exoRnase_RnaIDs': ['EG11620_RNA[c]', 'G7175_RNA[c]',
				'EG10858_RNA[c]', 'EG10863_RNA[c]', 'EG11259_RNA[c]',
				'EG11547_RNA[c]', 'EG10746_RNA[c]', 'G7842_RNA[c]',
				'EG10743_RNA[c]'],

			'rProteins': ['EG10872-MONOMER[c]', 'EG10879-MONOMER[c]',
				'EG11232-MONOMER[c]', 'EG10877-MONOMER[c]', 'EG10876-MONOMER[c]',
				'EG10892-MONOMER[c]', 'EG10874-MONOMER[c]',	'EG50001-MONOMER[c]',
				'EG10875-MONOMER[c]', 'EG10884-MONOMER[c]', 'EG11231-MONOMER[c]',
				'EG10887-MONOMER[c]', 'EG10871-MONOMER[c]',	'EG10878-MONOMER[c]',
				'EG10886-MONOMER[c]', 'EG10870-MONOMER[c]', 'EG10889-MONOMER[c]',
				'EG10891-MONOMER[c]', 'EG10888-MONOMER[c]', 'EG50002-MONOMER[c]',
				'EG10869-MONOMER[c]', 'EG10882-MONOMER[c]',	'EG10883-MONOMER[c]',
				'EG10885-MONOMER[c]', 'EG10890-MONOMER[c]',	'EG10864-MONOMER[c]',
				'EG10881-MONOMER[c]', 'EG10865-MONOMER[c]', 'EG10868-MONOMER[c]',
				'EG10880-MONOMER[c]', 'EG10867-MONOMER[c]', 'EG10873-MONOMER[c]',
				'EG10866-MONOMER[c]', 'EG10912-MONOMER[c]', 'EG10916-MONOMER[c]',
				'EG10920-MONOMER[c]', 'EG10914-MONOMER[c]', 'EG10909-MONOMER[c]',
				'EG10903-MONOMER[c]', 'EG10911-MONOMER[c]', 'EG10904-MONOMER[c]',
				'EG10900-MONOMER[c]', 'EG10901-MONOMER[c]', 'EG10905-MONOMER[c]',
				'EG10915-MONOMER[c]', 'EG10918-MONOMER[c]', 'EG10919-MONOMER[c]',
				'EG10907-MONOMER[c]', 'EG11508-MONOMER[c]', 'EG10908-MONOMER[c]',
				'EG10906-MONOMER[c]', 'EG10910-MONOMER[c]', 'EG10902-MONOMER[c]',
				'EG10917-MONOMER[c]', 'EG10913-MONOMER[c]'],
		}

		bulkMoleculesBinomialDivision = createIdsWithCompartments(raw_data.metabolites)

		bulkMoleculesBinomialDivision.extend(createIdsWithCompartments(raw_data.rnas))
		bulkMoleculesBinomialDivision.extend(createIdsWithCompartments(raw_data.proteins))
		bulkMoleculesBinomialDivision.extend(createIdsWithCompartments(raw_data.proteinComplexes))
		bulkMoleculesBinomialDivision.extend(createIdsWithCompartments(raw_data.modifiedForms))
		bulkMoleculesBinomialDivision.extend(createIdsWithCompartments(raw_data.water))

		moleculeGroups['bulkMoleculesBinomialDivision'] = bulkMoleculesBinomialDivision
		moleculeGroups['bulkMoleculesEqualDivision'] = []

		self.__dict__.update(moleculeGroups)
