import subprocess
import os
import argparse

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


def main(out_directory):
	allSimulations = findDirectories(out_directory)

	for sim_dir in allSimulations:
		all_seeds = findDirectories(sim_dir)
		all_seeds.pop(all_seeds.index(os.path.join(sim_dir, 'kb')))
		all_seeds.pop(all_seeds.index(os.path.join(sim_dir, 'metadata')))

		for seed_idx, seed in enumerate(all_seeds):
			plotDir = os.path.join(seed, 'plotOut')
			print 'Working on {}/{}: {}'.format(seed_idx,len(all_seeds),seed)
			if os.path.exists(plotDir):
				allPdfPlots = findFiles(plotDir, '.pdf')
				if allPdfPlots:
					if not os.path.exists(os.path.join(plotDir, 'pdf_plots')):
						os.mkdir(os.path.join(plotDir, 'pdf_plots'))
					if not os.path.exists(os.path.join(plotDir, 'low_res_plots')):
						os.mkdir(os.path.join(plotDir, 'low_res_plots'))

					for idx,pdfPlot in enumerate(allPdfPlots):
						plotName = pdfPlot[:-4]
						print '{}/{} Converting {} to {}'.format(idx, len(allPdfPlots), plotName + ".pdf", plotName + ".svg")
						subprocess.call(
								[
									"inkscape",
									"-l",
									os.path.join(plotDir, plotName + ".svg"),
									os.path.join(plotDir, plotName + ".pdf")
								]
							)

						print 'Creating low resolution {}'.format(plotName + ".png")
						subprocess.call(
								[
									"inkscape",
									"--export-png",
									os.path.join(plotDir, "low_res_plots", plotName + ".png"),
									"--export-dpi=30",
									os.path.join(plotDir, plotName + ".svg")
								]
							)

						print 'Moving {} to {}'.format(plotName + ".pdf", "pdf_plots")
						subprocess.call(
								[
									"mv",
									os.path.join(plotDir, plotName + ".pdf"),
									os.path.join(plotDir, "pdf_plots", plotName + ".pdf")
								]
							)


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("out_directory", help = "wcEcoli/out directory", type = str)

	args = parser.parse_args().__dict__

	main(args["out_directory"])