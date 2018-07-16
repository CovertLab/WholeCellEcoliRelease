"""Test the memory_debug utility.

Running it this way reveals the print output with Node IDs:
	python -m wholecell.tests.utils.test_memory_debug

In any case, you should see GC messages like:
	gc: uncollectable <Node 0x110118050>
	gc: uncollectable <Node 0x110118090>
	gc: uncollectable <Node 0x110118110>
	gc: uncollectable <dict 0x1101124b0>
	gc: uncollectable <dict 0x110112a28>

# The test case is adapted from https://pymotw.com/2/gc/
"""

from __future__ import absolute_import
from __future__ import division

import gc
import nose.plugins.attrib as noseAttrib
import unittest

from wholecell.utils import memory_debug


class Node(object):
	def __init__(self, name):
		self.name = str(name)
		self.link = None
		self.link2 = None

	def __repr__(self):
		return "{}{} 0x{:x} dict 0x{:x}".format(
			type(self).__name__, self.name, id(self), id(self.__dict__))

	def __del__(self):
		print "    {} __del__()".format(self)


class Test_memory_debug(unittest.TestCase):
	@noseAttrib.attr('smalltest')
	def test_memory_debug(self):
		precount = len(gc.garbage)

		with memory_debug.detect_leaks(enabled=True):
			nodes = [Node(i) for i in xrange(6)]

			# N0 -> N1 -> N2 are not in a cycle and should be collectable.
			nodes[0].link = nodes[1]
			nodes[1].link = nodes[2]

			# Make a cycle + spur N3 <-> N4 -> N5 of uncollectable Nodes.
			# N3's and N4's dicts will be in the cycle.
			nodes[3].link = nodes[4]
			nodes[4].link = nodes[3]
			nodes[4].link2 = nodes[5]

			print "Dropping: Should __del__ {}.".format(nodes[:3])
			uncollectable = str(nodes[3:])  # don't retain the Nodes
			nodes = []

			print ("Collecting: {} and some of their dicts should be"
				   " uncollectable.").format(uncollectable)
			# Why is Node5's dict collectable?

		# gc.garbage holds Node3 .. Node5.
		self.assertEqual(precount + 3, len(gc.garbage),
			'Uncollectable: {}'.format(gc.garbage))


if __name__ == '__main__':
	unittest.main()
