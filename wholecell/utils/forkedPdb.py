from __future__ import absolute_import, division, print_function

import sys
from IPython.core.debugger import Pdb

class ForkedPdb(Pdb):
	""" Usage:
		==========
		from wholecell.utils import ForkedPdb
		ForkedPdb().set_trace()
		==========
		Don't forget to actually instantiate an
		instance of ForkedPdb when calling set_trace()

		Note: Python 3 has a global function breakpoint().
	"""
	def set_trace(self, frame = None):
		_stdin = sys.stdin
		sys.stdin = open("/dev/stdin")
		if frame is None:
			# noinspection PyUnresolvedReferences
			frame = sys._getframe().f_back
		Pdb("Linux").set_trace(frame)
