import unittest
import subprocess
import json

import os

TEMP_FILE_NAME = 'fixture_tmp.json'
import wholecell.utils.config
N_PROCESSES = wholecell.utils.config.N_PROCESSES

class BaseLargeTest(unittest.TestCase):
	fixtureOpts = None # assigned by subclass

	@classmethod
	def setUpClass(cls):
		generateFixtures = True
		if generateFixtures:
			# Check that fixtureOpts is populated and write to temporary file
			if cls.fixtureOpts is None:
				raise Exception('No fixture options defined for test which inherits ' + str(self) + '!')

			json.dump(cls.fixtureOpts, open(TEMP_FILE_NAME, 'w'))
			
			# Load json file into generateLargeTestFixtures
			# Use info to run simulation with MPI
			returnCode = subprocess.call(
				['mpirun', '-n', str(N_PROCESSES), 'python', 'wholecell/tests/processes/generate_large_test_fixtures.py']
				)

			if returnCode != 0:
				# TODO: confirm that nonzero return codes are errors
				raise Exception('mpirun call failed with return code {}'.format(returnCode))
			
			# Delete json file
			os.remove(TEMP_FILE_NAME)
