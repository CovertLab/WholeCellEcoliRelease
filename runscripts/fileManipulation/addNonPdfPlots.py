from __future__ import absolute_import, division, print_function

import os
import argparse

from wholecell.utils import constants
import wholecell.utils.filepath as fp


def findDirectories(directory):
	allFiles = os.listdir(directory)
	onlyDirs = []
	for f in allFiles:
		temp = os.path.join(directory,f)
		if os.path.isdir(temp): onlyDirs.append(temp) 
	return onlyDirs

def findFiles(directory,typeFile):
	if not os.path.isdir(directory): return []
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
		all_seeds.pop(all_seeds.index(os.path.join(sim_dir, constants.KB_DIR)))
		all_seeds.pop(all_seeds.index(os.path.join(sim_dir, constants.METADATA_DIR)))

		for seed_idx, seed in enumerate(all_seeds):
			plotDir = os.path.join(seed, constants.PLOTOUT_DIR)
			print('Working on {}/{}: {}'.format(seed_idx,len(all_seeds),seed))
			if os.path.exists(plotDir):
				allPdfPlots = findFiles(plotDir, '.pdf')
				if allPdfPlots:
					if not os.path.exists(os.path.join(plotDir, 'pdf_plots')):
						os.mkdir(os.path.join(plotDir, 'pdf_plots'))
					if not os.path.exists(os.path.join(plotDir, 'low_res_plots')):
						os.mkdir(os.path.join(plotDir, 'low_res_plots'))

					for idx,pdfPlot in enumerate(allPdfPlots):
						plotName = pdfPlot[:-4]
						print('{}/{} Converting {} to {}'.format(idx, len(allPdfPlots), plotName + ".pdf", plotName + ".svg"))
						fp.run_cmd(
								[
									"inkscape",
									"-l",
									os.path.join(plotDir, plotName + ".svg"),
									os.path.join(plotDir, plotName + ".pdf")
								]
							)

						print('Creating low resolution {}'.format(plotName + ".png"))
						fp.run_cmd(
								[
									"inkscape",
									"--export-png",
									os.path.join(plotDir, "low_res_plots", plotName + ".png"),
									"--export-dpi=30",
									os.path.join(plotDir, plotName + ".svg")
								]
							)

						print('Moving {} to {}'.format(plotName + ".pdf", "pdf_plots"))
						fp.run_cmd(
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
