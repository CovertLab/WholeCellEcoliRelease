import argparse
import os
import subprocess

def findDirectories(directory):
	allFiles = os.listdir(directory)
	onlyDirs = []
	for f in allFiles:	
		temp = os.path.join(directory,f)
		if os.path.isdir(temp): onlyDirs.append(temp) 
	return onlyDirs

def findFiles(directory,typeFile):
	if os.path.isdir(directory) == False: return []
	allFiles = os.listdir(directory)
	onlyFiles = []
	for f in allFiles:	
		temp = os.path.join(directory,f)
		if os.path.isfile(temp):
			if typeFile in f:
				onlyFiles.append(f) 
	return onlyFiles

def main(out_dir):
	allSimulations = findDirectories(out_dir)

	for simPath in allSimulations:

		# Checks for old simulation format with 000000 seed
		if not os.path.exists(os.path.join(simPath, "000000")):
			return

		# Create list of folders to move down one level
		subFolders = findDirectories(simPath)
		metadataPath = os.path.join(simPath, "metadata")
		kbPath = os.path.join(simPath, "kb")
		subFolders.pop(subFolders.index(metadataPath))
		subFolders.pop(subFolders.index(kbPath))

		# Create variant level
		subprocess.call([
			'mkdir',
			os.path.join(simPath, 'wildtype_000000')
			])
		subprocess.call([
			'mkdir',
			os.path.join(simPath, 'wildtype_000000', 'kb')
			])
		subprocess.call([
			'mkdir',
			os.path.join(simPath, 'wildtype_000000', 'metadata')
			])

		with open(os.path.join(simPath, 'wildtype_000000', 'metadata', 'description'), 'w') as writefile:
			writefile.write('Control simulation')

		with open(os.path.join(simPath, 'wildtype_000000', 'metadata', 'short_name'), 'w') as writefile:
			writefile.write('control')

		# Move files down to variant level
		for subFolder in subFolders:
			subprocess.call([
				'mv',
				subFolder,
				os.path.join(simPath, 'wildtype_000000')
				])

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("out_directory", help = "wcEcoli/out directory", type = str)

	args = parser.parse_args().__dict__

	main(args["out_directory"])