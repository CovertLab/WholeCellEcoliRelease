from __future__ import annotations

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"doubling_time_histogram.py",
	"excess_protein_monomers.py",
	"growth_histograms.py",
	"mRNA_copy_numbers.py",
	"massFractionComparison.py",
	"mRNA_copy_numbers.py",
	"mRNA_length_histogram.py",
	"polycistronic_transcription.py",
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"growth_histograms.py",
		"massFractionComparison.py",
		],
	}
