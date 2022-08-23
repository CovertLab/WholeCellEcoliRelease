"""
Knockout expression of mechanistic amino acid synthesis genes in rich conditions.

Modifies:
	attributes from condition variant
	attributes from gene_knockout variant

Expected variant indices (depends on the enzymes in sim_data.process.metabolism.synthesis_enzymes):
	0: control
	1-31: enzyme to knockout
"""

from .condition import condition


def aa_synthesis_ko(sim_data, index):
	_, sim_data = condition(sim_data, 2)

	metabolism = sim_data.process.metabolism
	complexation = sim_data.process.complexation
	monomer_data = sim_data.process.translation.monomer_data
	cistron_data = sim_data.process.transcription.cistron_data
	gene_data = sim_data.process.replication.gene_data

	# Enzymes involved in mechanistic amino acid synthesis
	synthesis_enzymes = metabolism.aa_enzymes[metabolism.enzyme_to_amino_acid_fwd.sum(1).astype(bool)]
	synthesis_monomers = sorted({
		subunit
		for enzyme in synthesis_enzymes
		for subunit in complexation.get_monomers(enzyme)['subunitIds']
		})

	# Map monomers to RNA for a knockout
	monomer_to_cistron = {monomer['id']: monomer['cistron_id'] for monomer in monomer_data}
	cistron_to_index = {cistron['id']: i for i, cistron in enumerate(cistron_data)}
	cistron_to_symbol = {gene['cistron_id']: gene['symbol'] for gene in gene_data}
	cistrons = [
		monomer_to_cistron[monomer]
		for monomer in synthesis_monomers
		if monomer in monomer_to_cistron
		]

	n_variants = len(cistrons) + 1
	if index > n_variants:
		raise ValueError(f'Variant index {index} is not supported. Choose between 0 and {n_variants}')

	if index > 0:
		cistron = cistrons[index - 1]
		rna_index = cistron_to_index[cistron]
		sim_data.adjust_final_expression([rna_index], [0])

		symbol = cistron_to_symbol[cistron]
		name=f'{symbol} KO'
		desc=f'{symbol} KO in rich media'
	else:
		name='control'
		desc='Control in rich media'

	return dict(shortName=name, desc=desc), sim_data
