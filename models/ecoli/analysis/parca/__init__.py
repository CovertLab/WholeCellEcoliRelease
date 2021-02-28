from __future__ import absolute_import, division, print_function

# Active analysis modules to run, in this order.
# Tip: Edit this during development to run ones you're working on first.
ACTIVE = [
	'expression_probabilities.py',
	'fold_changes.py',
	'interpolation.py',
	'metabolite_concentrations.py',
	'ppgpp_expression.py',
	'tf_target.py',
]

TAGS = {
	'ACTIVE': ACTIVE,   # all active analyses in this category
	'CORE': [           # the default list to run in development
		'metabolite_concentrations.py',
		'interpolation.py',
		],
	}
