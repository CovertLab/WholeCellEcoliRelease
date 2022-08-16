from __future__ import annotations

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"coexpression_probabilities.py",
	"doubling_time_histogram.py",
	"excess_protein_monomers.py",
	"gene_position_vs_expression_change.py",
	"growth_histograms.py",
	"massFractionComparison.py",
	"mRNA_copy_numbers.py",
	"mRNA_copy_numbers_growth_genes.py",
	"mRNA_counts_histogram.py",
	"mRNA_length_histogram.py",
	"mRNA_mass_histogram.py",
	"polycistronic_transcription.py",
	"protein_mass_histogram.py",
	"protein_stoichiometry.py",
	"proteomics_fluxomics_comparison.py",
	"proteomics_fluxomics_validation.py",
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"growth_histograms.py",
		"massFractionComparison.py",
		],
	}
