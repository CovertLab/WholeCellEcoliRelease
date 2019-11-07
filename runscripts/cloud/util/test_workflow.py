#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

from runscripts.cloud.util.workflow_cli import WorkflowCLI


class TestWorkflow(WorkflowCLI):
	"""A test workflow for integration and regression tests of the cloud
	workflow software.

	The run will never reach "WORKFLOW COMPLETE" since some tasks intentionally
	fail and another task will never run for lack of inputs from failed tasks.
	"""

	DEFAULT_TIMEOUT = 30  # in seconds

	def build(self, args):
		"""Build the workflow."""
		lines_filename = '/tmp/lines.txt'
		code = (
			"with open('" + lines_filename + "', 'w') as f:\n"
			"  for i in range(10):\n"
			"    f.write('This is line {}\\n'.format(i))\n"
			"    print('line {}'.format(i))")
		self.add_task(
			name='lines',
			outputs=[lines_filename],
			command=['python', '-u', '-c', code])

		self.add_task(
			name='count',
			inputs=[lines_filename],
			command=['wc', lines_filename])

		# Expected:  "wc: /tmp/lines.txt: No such file or directory"
		# because this task spec didn't request the input file.
		error_no_such_file_out = '/tmp/expected_no_such_file.txt'
		self.add_task(
			name='expected_no_such_file',
			inputs=[],
			outputs=['>' + error_no_such_file_out],
			command=['wc', lines_filename])

		# Expected:  "IndexError: tuple index out of range"
		index_error_out = '/tmp/expected_index_out_of_range.txt'
		self.add_task(
			name='expected_index_out_of_range',
			outputs=['>' + index_error_out],
			command=['python', '-u', '-c', "()[1]"])

		# This task depends on error task outputs as a regression test to check
		# that the workflow engine doesn't keep retrying any of them.
		# Expected:  This task never runs since its inputs never arrive.
		self.add_task(
			name='expected_to_never_run',
			inputs=[error_no_such_file_out, index_error_out],
			command=['cat', error_no_such_file_out, index_error_out])

		# This task writes files into an output dir to test the file ownership
		# of files created by the process inside the Docker container, not by
		# the Sisyphus worker (which creates files and directories explicitly
		# named in task `inputs` and `outputs`).
		# Expected:  The text files on the worker server have ordinary user and
		# group ownership, not root, so the worker can delete them without error.
		output_dir = '/tmp/output/dir/'
		code = (
			"for i in range(4):\n"
			"  name = '" + output_dir + "{}.txt'.format(i)\n"
			"  print('Wrote ' + name)\n"
			"  with open(name, 'w') as f:\n"
			"    f.write('This is file {}\\n'.format(i))\n")
		self.add_task(
			name='output_dir',
			inputs=(),
			outputs=(output_dir,),
			command=['python', '-u', '-c', code])

		# Download and append to a file written by a previous Task to test file
		# permissions, e.g. it's not owned by root so the task can overwrite it.
		# It wouldn't fit Gaia's functional data flow model to output a file
		# back to storage with an input's filename since that means ambiguous
		# responsibility for which task creates it. This test skirts that by
		# printing the appended file and uploading stdout to prove that it
		# succeeded.
		code = (
			"fn = '" + output_dir + "1.txt'\n"
			"with open(fn, 'a') as f:\n"
			"  f.write('This is still file 1\\n')\n"
			"with open(fn, 'r') as f:\n"
			"  print(f.read())\n")
		self.add_task(
			name='overwrite',
			inputs=(output_dir,),
			outputs=['>/tmp/overwrite.txt'],
			command=['python', '-u', '-c', code])

		# test a timeout
		code = (
			"from time import sleep\n"
			"for i in range(100):\n"
			"  sleep(1)\n"
			"  print('{:3} seconds'.format(i))")
		self.add_task(
			name='expected_timeout',
			timeout=8,
			command=['python', '-u', '-c', code])


if __name__ == '__main__':
	TestWorkflow().cli()
