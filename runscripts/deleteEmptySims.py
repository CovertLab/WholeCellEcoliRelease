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
	print 'ONLY WORKS ON OLD FORMAT SIMULATIONS!'
	allSimulations = findDirectories(out_dir)
	to_delete = []
	for sim_idx, sim_dir in enumerate(allSimulations):
		markDelete = False
		print 'Working on {}'.format(sim_dir)
		
		all_seeds = findDirectories(sim_dir)
		if len(all_seeds) == 0:
			# Delete if no simulations were run
			markDelete = True
		all_seeds.pop(all_seeds.index(os.path.join(sim_dir, 'kb')))
		all_seeds.pop(all_seeds.index(os.path.join(sim_dir, 'metadata')))
		
		if not markDelete:
			for seed in all_seeds:
				plotDir = os.path.join(seed, 'plotOut')
				simDir = os.path.join(seed, 'simOut')
				if not os.path.exists(plotDir):
					# Delete if no plot dir
					markDelete = True
				if not os.path.exists(simDir):
					# Delete if no sim dir
					markDelete = True
			
		if markDelete:
			to_delete.append(sim_idx)
	print 'ONLY WORKS ON OLD FORMAT SIMULATIONS!'

	for idx in to_delete:
		subprocess.call([
			"rm",
			"-fr",
			allSimulations[idx]
		])		

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("out_directory", help = "wcEcoli/out directory", type = str)

	args = parser.parse_args().__dict__

	main(args["out_directory"])