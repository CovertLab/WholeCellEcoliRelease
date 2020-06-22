from __future__ import absolute_import, division, print_function

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	'interpolation.py',
	'metabolite_concentrations.py',
	'tf_target.py',
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		'metabolite_concentrations.py',
		'interpolation.py',
		],
	}
