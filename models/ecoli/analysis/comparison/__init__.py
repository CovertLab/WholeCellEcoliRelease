from __future__ import annotations

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	"growth_histograms.py",
	"massFractionComparison.py",
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		"growth_histograms.py",
		"massFractionComparison.py",
		],
	}
