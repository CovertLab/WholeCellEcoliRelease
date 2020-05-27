'''

See README.md for details.

'''

from __future__ import absolute_import, division, print_function

import os

from .model import simulate

OUTDIR = os.path.join(
	os.path.split(__file__)[0],
	'out'
	)

DEMONSTRATIONS = (
	dict(),
	dict(
		do_longstep = True
		),
	dict(
		do_integral = True
		),
	dict(
		do_integral = True,
		do_bootstrap = True
		),
	dict(
		do_longstep = True,
		do_integral = True,
		do_bootstrap = True
		),
	dict(
		do_longstep = True,
		do_integral = True,
		do_bootstrap = True,
		do_cut = True
		),
	dict(
		do_longstep = True,
		do_integral = True,
		do_bootstrap = True,
		do_cut = True,
		do_saturate = True
		)
	)

if not os.path.exists(OUTDIR):
	os.mkdir(OUTDIR)

for i, kwargs in enumerate(DEMONSTRATIONS):
	simulate(
		save_as = os.path.join(OUTDIR, 'demo{}.png'.format(i+1)),
		**kwargs
		)
