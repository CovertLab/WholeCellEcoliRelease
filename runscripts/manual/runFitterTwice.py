import os
import cPickle

from runscripts.manual.runFitter import RunFitter
from runscripts.reflect.object_tree import object_tree, diff_trees

if __name__ == '__main__':
	sim_data_file = os.path.join(
		os.getcwd(),
		'out',
		'manual',
		'kb',
		'simData_Most_Fit.cPickle')

	first = RunFitter()
	first.cli()
	with open(sim_data_file, "rb") as f:
		once = cPickle.load(f)
		
	second = RunFitter()
	second.cli()
	with open(sim_data_file, "rb") as f:
		twice = cPickle.load(f)

	a = object_tree(once)
	b = object_tree(twice)
	aminusb = diff_trees(a, b)
	bminusa = diff_trees(b, a)

	print(aminusb)
	print(bminusa)
