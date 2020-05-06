from __future__ import absolute_import, division, print_function

from vivarium.environment.boot import BootEnvironment
from vivarium.environment.boot import wrap_boot, wrap_init_basic
from environment.wcecoli_process import wcEcoliAgent

class BootEcoli(BootEnvironment):
	def __init__(self):
		super(BootEcoli, self).__init__()
		self.agent_types.update({
			'ecoli': wrap_boot(wrap_init_basic(wcEcoliAgent), {
                'volume': 1.0})})

def run():
	boot = BootEcoli()
	boot.execute()

if __name__ == '__main__':
	run()
