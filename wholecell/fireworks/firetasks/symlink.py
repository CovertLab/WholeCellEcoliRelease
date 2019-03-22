import os
import time

from fireworks import FireTaskBase, explicit_serialize


@explicit_serialize
class SymlinkTask(FireTaskBase):

	_fw_name = "SymlinkTask"
	required_params = ["to", "link"]
	optional_params = ["overwrite_if_exists"]

	def run_task(self, fw_spec):
		print "%s: Creating symlink" % (time.ctime())

		if self["overwrite_if_exists"]:
			if os.path.exists(self["link"]):
				os.unlink(self["link"])

		os.symlink(self["to"], self["link"])
