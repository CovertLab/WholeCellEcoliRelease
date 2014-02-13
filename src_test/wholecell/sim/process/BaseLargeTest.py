import unittest
import subprocess
import json

import os

TEMP_FILE_NAME = 'fixture_tmp.json'

class BaseLargeTest(unittest.TestCase):
	@classmethod
	def setUpClass(cls):
		generateFixtures = True
		if generateFixtures:
			# Check that fixture_opts is populated and write to temporary file
			if not len(cls.fixture_opts):
				raise Exception, 'No fixture options defined for test which inherits ' + str(self) + '!\n'
			else:
				outfile = open(TEMP_FILE_NAME, 'w')
				outfile.write(json.dumps(cls.fixture_opts))
				outfile.close()
			
			# TODO: Load n from config file
			n = 4
			# Load json file into generateLargeTestFixtures
			# Use info to run simulation with MPI
			subprocess.call(['mpirun', '-n', str(n), 'python', 'src_test/wholecell/sim/process/generateLargeTestFixtures.py'])

			# Delete json file
			os.remove(TEMP_FILE_NAME)