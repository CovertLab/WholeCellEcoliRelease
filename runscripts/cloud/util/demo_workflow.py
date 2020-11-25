#!/usr/bin/env python

from runscripts.cloud.util.workflow_cli import WorkflowCLI


class DemoWorkflow(WorkflowCLI):
	"""A demo workflow for quick demonstration of how to build a workflow."""

	DEFAULT_TIMEOUT = 10  # in seconds

	def build(self, args):
		"""Build the workflow."""
		lines_filename = self.internal('lines.txt')
		code = ("with open('" + lines_filename + "', 'w') as f:\n"
			"  for i in range(100):\n"
			"    f.write('This is line {}\\n'.format(i))\n"
			"    print('hello {}'.format(i))")
		self.add_task(
			name='lines',
			outputs=[lines_filename],
			command=['python', '-u', '-c', code])

		self.add_task(
			name='count',
			inputs=[lines_filename],
			outputs=['>' + self.internal('count.txt')],
			command=['wc', lines_filename])


if __name__ == '__main__':
	DemoWorkflow().cli()
