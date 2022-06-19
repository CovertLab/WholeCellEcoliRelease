"""
Knockout expression of mechanistic amino acid synthesis genes and then shift from rich to minimal.
Selected genes are not essential in minimal media because of redundancy.

Modifies:
	attributes from remove_aas_shift variant
	attributes from adjust_final_expression

Expected variant indices (depends on SYNTHESIS_GENES):
	0: control shift from rich to minimal without a knockout
	1-7: enzyme knockouts shifted to minimal
	8-14: enzyme knockouts shifted to minimal plus corresponding AA
"""

from .remove_aas_shift import remove_aas_shift


# TODO: generalize index based on remove_aas_shift file
SYNTHESIS_GENES = {
	'alaC': 4, 'alaA': 4,  # Ala synthesis
	'asnA': 6, 'asnB': 6,  # Asn synthesis
	'gltB': 8, 'gltD': 8, 'gdhA': 8,  # Glt synthesis
	}
SHIFT_TIME = 2 * 3600  # 2 hrs


def aa_synthesis_ko_shift(sim_data, index):
	# Map genes to RNA for a knockout
	replication = sim_data.process.replication
	transcription = sim_data.process.transcription
	symbol_to_cistron = {gene['symbol']: gene['cistron_id'] for gene in replication.gene_data}
	cistron_to_index = {cistron['id']: i for i, cistron in enumerate(transcription.cistron_data)}

	n_genes = len(SYNTHESIS_GENES)
	n_variants = 2 * n_genes + 1
	if index >= n_variants:
		raise ValueError(f'Variant index {index} is not supported. Choose between 0 and {n_variants-1}')

	gene = list(SYNTHESIS_GENES.keys())[(index - 1) % n_genes]
	shift_index = SYNTHESIS_GENES[gene] if index > n_genes else 23  # TODO: not hardcoded index

	# Shift to new media
	# TODO: make this less hacky (make a function that is used by both remove_aas_shift and here
	# to make a shift and set a shift time)
	old_timelines = set(sim_data.external_state.saved_timelines.keys())
	shift_desc, sim_data = remove_aas_shift(sim_data, shift_index)
	new_timeline_id = list(set(sim_data.external_state.saved_timelines.keys()) - old_timelines)[0]
	new_timeline = sim_data.external_state.saved_timelines[new_timeline_id]
	new_timeline[1] = (SHIFT_TIME, new_timeline[1][1])

	if index > 0:
		rna_index = cistron_to_index[symbol_to_cistron[gene]]
		sim_data.adjust_final_expression([rna_index], [0])  # 0 for KO

		name=f'{shift_desc["shortName"]}-{gene} KO'
		desc=f'{shift_desc["desc"]} with {gene} KO'
	else:
		name='control'
		desc='Control shift from rich to minimal'

	return dict(shortName=name, desc=desc), sim_data
