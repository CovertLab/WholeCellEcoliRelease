#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

from runscripts.cloud.util.workflow_cli import WorkflowCLI


class TestWorkflow(WorkflowCLI):
	"""A test workflow for integration and regression tests of the workflow software."""

	def build(self, args):
		# Build the workflow.
		lines_filename = '/tmp/lines.txt'
		code = ("with open('" + lines_filename + "', 'w') as f:\n"
			"  for i in range(10):\n"
			"    f.write('This is line {}\\n'.format(i))\n"
			"    print 'hello', i")
		self.add_task(
			name='lines',
			outputs=[lines_filename],
			command=['python', '-u', '-c', code])

		self.add_task(
			name='count',
			inputs=[lines_filename],
			outputs=['>/tmp/count.log'],
			command=['wc', lines_filename])

		# Expected:  wc: /tmp/lines.txt: No such file or directory
		# because we "forgot" to download the input file.
		error_test_log = '/tmp/error-test.log'
		self.add_task(
			name='error_test',
			inputs=[],
			outputs=['>' + error_test_log],
			command=['wc', lines_filename])

		# Expected:  IndexError: tuple index out of range
		exception_log = '/tmp/index_exception.log'
		self.add_task(
			name='index_exception',
			outputs=['>' + exception_log],
			command=['python', '-u', '-c', "()[1]"])

		# This task depends on error tasks to test that they don't keep retrying.
		# Expected:  Normal output concatenating the two input logs with exceptions.
		# This is a regression test.
		two_logs = '/tmp/two.log'
		self.add_task(
			name='error_watcher',
			inputs=[error_test_log, exception_log],
			outputs=['>' + two_logs],
			command=['cat', error_test_log, exception_log])

		# This task writes files into an output dir to test the file ownership
		# of files created inside the Docker container, not by the Sisyphus
		# worker (which creates named input and output files and directories).
		# Expected:  The text files on the worker server have ordinary user and
		# group ownership, not root, so the worker can delete them without error.
		output_dir = '/tmp/output/dir/'
		code = ("for i in range(4):\n"
			"  name = '" + output_dir + "{}.txt'.format(i)\n"
			"  print name\n"
			"  with open(name, 'w') as f:\n"
			"    f.write('This is file {}\\n'.format(i))\n")
		self.add_task(
			name='output_dir',
			inputs=(),
			outputs=(output_dir,),
			command=['python', '-u', '-c', code])

		# Download and overwrite a file written by a previous Task to test file
		# permissions, e.g. not owned by root. It wouldn't fit Gaia's functional
		# data flow model to upload a file back to Gaia since that means
		# ambiguous responsibility for creating it. This test skirts that by
		# uploading stdout, not the overwritten file.
		code = ("fn = '" + output_dir + "1.txt'\n"
			"with open(fn, 'a') as f:\n"
			"  f.write('This is still file 1\\n')\n"
			"with open(fn, 'r') as f:\n"
			"  print(f.read())\n")
		self.add_task(
			name='overwrite',
			inputs=(output_dir,),
			outputs=['>/tmp/overwrite.log'],
			command=['python', '-u', '-c', code])


if __name__ == '__main__':
	TestWorkflow().cli()
