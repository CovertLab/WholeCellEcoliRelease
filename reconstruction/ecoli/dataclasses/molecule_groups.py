"""
SimulationData moleculeGroups
"""

from __future__ import absolute_import, division, print_function

POLYMERIZED_FRAGMENT_PREFIX = 'polymerized_'


class MoleculeGroups(object):
	"""
	Helper class to extract molecule IDs of "special" groups of molecules. All
	values returned are lists of strings.
	"""

	def __init__(self, raw_data, sim_data):
		self._build_molecule_groups(sim_data)

	def _build_molecule_groups(self, sim_data):
		aa_ids = list(sim_data.amino_acid_code_to_id_ordered.values())
		ntp_ids = list(sim_data.ntp_code_to_id_ordered.values())
		dntp_ids = list(sim_data.dntp_code_to_id_ordered.values())
		polymerized_aa_ids = [
			POLYMERIZED_FRAGMENT_PREFIX + aa_id for aa_id in aa_ids]
		polymerized_ntp_ids = [
			POLYMERIZED_FRAGMENT_PREFIX + ntp_id for ntp_id in ntp_ids]
		polymerized_dntp_ids = [
			POLYMERIZED_FRAGMENT_PREFIX + dntp_id for dntp_id in dntp_ids]

		codon_ids = []
		ntp_abbreviations = [ntp[0] for ntp in ntp_ids]
		for nucleotide_0 in ntp_abbreviations:
			for nucleotide_1 in ntp_abbreviations:
				for nucleotide_2 in ntp_abbreviations:
					codon = nucleotide_0 + nucleotide_1 + nucleotide_2

					# Skip stop codons (except UGA, which codes for selenocysteine)
					if codon in ['UAA', 'UAG']:
						continue

					codon_ids.append(codon)

		molecule_groups = {
			'amino_acids': aa_ids,
			'ntps': ntp_ids,
			'dntps': dntp_ids,

			'polymerized_amino_acids': polymerized_aa_ids,
			'polymerized_ntps': polymerized_ntp_ids,
			'polymerized_dntps': polymerized_dntp_ids,
			'polymerized_subunits': polymerized_aa_ids + polymerized_ntp_ids + polymerized_dntp_ids,

			'codons': codon_ids,

			'initiator_trnas': ['RNA0-306[c]', 'metY-tRNA[c]', 'metZ-tRNA[c]', 'metW-tRNA[c]'],
			'elongator_trnas': ['metT-tRNA[c]', 'metU-tRNA[c]'],

			# EG11508-MONOMER[c] (sra) is removed from the list because it is
			# not essential and should have very low counts in rich media
			's30_proteins':	['EG10912-MONOMER[c]', 'EG10916-MONOMER[c]',
				'EG10906-MONOMER[c]', 'EG10914-MONOMER[c]', 'EG10909-MONOMER[c]',
				'EG10903-MONOMER[c]', 'EG10911-MONOMER[c]', 'EG10904-MONOMER[c]',
				'EG10900-MONOMER[c]', 'EG10901-MONOMER[c]', 'EG10905-MONOMER[c]',
				'EG10915-MONOMER[c]', 'EG10918-MONOMER[c]', 'EG10919-MONOMER[c]',
				'EG10907-MONOMER[c]', 'EG10908-MONOMER[c]', 'EG10920-MONOMER[c]',
				'EG10910-MONOMER[c]', 'EG10902-MONOMER[c]', 'EG10917-MONOMER[c]',
				'EG10913-MONOMER[c]'],
			's30_16s_rRNA': ['RRSA-RRNA[c]', 'RRSB-RRNA[c]', 'RRSC-RRNA[c]',
				'RRSD-RRNA[c]', 'RRSE-RRNA[c]', 'RRSG-RRNA[c]', 'RRSH-RRNA[c]'],

			's50_protein_complexes': ['CPLX0-3956[c]'],
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
			's50_23s_rRNA': ['RRLA-RRNA[c]', 'RRLB-RRNA[c]', 'RRLC-RRNA[c]',
				'RRLD-RRNA[c]', 'RRLE-RRNA[c]', 'RRLG-RRNA[c]', 'RRLH-RRNA[c]'],
			's50_5s_rRNA': ['RRFA-RRNA[c]', 'RRFB-RRNA[c]', 'RRFC-RRNA[c]',
				'RRFD-RRNA[c]', 'RRFE-RRNA[c]', 'RRFG-RRNA[c]', 'RRFH-RRNA[c]'],

			'lipids': ['CPD-8260[c]', 'CPD-12819[c]', 'CPD-12824[c]'],
			'polyamines': ['GAMMA-GLUTAMYL-PUTRESCINE[c]', 'PUTRESCINE[c]',
				'GLUTATHIONYLSPERMIDINE[c]', 'SPERMIDINE[c]',
				'N1-ACETYLSPERMINE[c]', 'SPERMINE[c]'],

			# TODO: 'EG10245-MONOMER[c]' (DNAP III subunit tau) should be added
			# 	to the list of trimer subunits once frame-shifting proteins are
			# 	produced.
			'replisome_trimer_subunits': ['CPLX0-2361[c]', 'CPLX0-3761[c]'],
			'replisome_monomer_subunits': ['CPLX0-3621[c]', 'EG10239-MONOMER[c]',
				'EG11500-MONOMER[c]', 'EG11412-MONOMER[c]'],

			'exoRNases': ['EG11620-MONOMER[c]', 'G7175-MONOMER[c]',
				'EG10858-MONOMER[c]', 'EG10863-MONOMER[c]', 'EG11259-MONOMER[c]',
				'EG11547-MONOMER[c]', 'EG10746-MONOMER[c]', 'G7842-MONOMER[c]',
				'EG10743-MONOMER[c]'],
			'endoRNase_rnas': ['EG10856_RNA', 'EG10857_RNA',
				'EG10859_RNA', 'EG10860_RNA', 'EG10861_RNA',
				'EG10862_RNA', 'EG11299_RNA', 'G7175_RNA',
				'G7365_RNA'],
			'exoRNase_rnas': ['EG11620_RNA', 'G7175_RNA',
				'EG10858_RNA', 'EG10863_RNA', 'EG11259_RNA',
				'EG11547_RNA', 'EG10746_RNA', 'G7842_RNA',
				'EG10743_RNA'],
			'RNAP_subunits': ['RPOB-MONOMER[c]', 'RPOC-MONOMER[c]',
				'EG10893-MONOMER[c]'],

			# EG11508-MONOMER[c] (sra) is removed from the list because it is
			# not essential and should have very low counts in rich media
			'ribosomal_proteins': ['EG10872-MONOMER[c]', 'EG10879-MONOMER[c]',
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
				'EG10907-MONOMER[c]', 'EG10908-MONOMER[c]', 'EG10906-MONOMER[c]',
				'EG10910-MONOMER[c]', 'EG10902-MONOMER[c]', 'EG10917-MONOMER[c]',
				'EG10913-MONOMER[c]'],

			'carbon_sources': ['GLC[p]', 'ACET[p]', 'SUC[p]'],
		}

		# Initialize molecule groups for how molecules are split between two
		# daughter cells at cell division (populated later by InternalState)
		molecule_groups['bulk_molecules_binomial_division'] = []
		molecule_groups['bulk_molecules_equal_division'] = []

		molecule_groups['unique_molecules_active_ribosome_division'] = []
		molecule_groups['unique_molecules_RNA_division'] = []
		molecule_groups['unique_molecules_domain_index_division'] = []
		molecule_groups['unique_molecules_chromosomal_segment_division'] = []

		self.__dict__.update(molecule_groups)
