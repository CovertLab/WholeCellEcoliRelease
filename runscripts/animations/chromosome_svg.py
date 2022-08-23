'''
Creates series of svg files that can be stitched together using avconv to create an animation of chromosome replication
Input: repAnimation.tsv file in same file location with 4 columns (time, number of chromosomes, fraction complete for first pair of replication forks,
fraction complete for second set of replication forks)
Output: time series of svg images to folder svg
'''

from __future__ import absolute_import, division, print_function

import io
import os

import numpy as np

from wholecell.io import tsv
from wholecell.utils import filepath as fp


def boilerplate_start(h):
#	h.write("<html>\n")
#	h.write("<body>\n")
	h.write("<svg width=\"100%\" height=\"100%\"  viewBox=\"0 0 210 825\">\n")

def boilerplate_end(h):
	h.write("</svg>\n")
#	h.write("</body>\n")
#	h.write("</html>")

def svg_draw(h, count, frac1, frac2):
	if frac1 < 0 or frac1 >= 1:
		raise Exception("frac1 is out of bounds")
	if frac2 < 0 or frac2 >= 0.5:
		raise Exception("frac2 is out of bounds")
	if frac2 > frac1:
		raise Exception("frac1 must be greater than frac2")
	if count <= 0 or count > 2:
		raise Exception("count must be 1 or 2")
	# xVal = 0.
	yVal1 = 0.
	useOuter = 0
	spaceBetween = 15.

	h.write(" <defs>\n")
	h.write("  <symbol id=\"chromosome\">\n")
	h.write("   <g transform=\"translate(105,105) scale(1,-1)\">\n")
	h.write("    <g style=\"stroke: black; stroke-width: 3\">\n")
	h.write("     <path d=\"M 100 0\n")
	h.write("              A 100 100, 0, 0, 1, -100 0\n")
	h.write("              A 100 100, 0, 0, 1, 100 0\n")
	h.write("             \" fill-opacity=\"0.0\"/>\n")
	if frac1 > 0:
		xVal = 100 * np.sin(frac1 * np.pi)
		yVal1 = 100 * np.cos(frac1 * np.pi)
		if frac1 > 0.5:
			useOuter = 1
		h.write("     <path d=\"M %f %f\n" % (xVal, yVal1))
		h.write("              A 100 100, 0, %d, 0, -%f %f\" fill-opacity=\"0.0\"/>\n" % (useOuter, xVal, yVal1))
	if frac2 > 0:
		c1Center = 200 * np.cos(frac1 * np.pi)
		c2Center = 200 * np.cos(frac2 * np.pi)
		xVal = 100 * np.sin(frac2 * np.pi)
		yVal = -100 * np.cos(frac2 * np.pi) + c1Center
		h.write("     <path d=\"M %f %f\n" % (xVal, yVal))
		h.write("              A 100 100, 0, 0, 1, -%f %f\" fill-opacity=\"0.0\"/>\n" % (xVal, yVal))
		yVal = -100 * np.cos(frac2 * np.pi) + c2Center
		h.write("     <path d=\"M %f %f\n" % (xVal, yVal))
		h.write("              A 100 100, 0, 0, 0, -%f %f\" fill-opacity=\"0.0\"/>\n" % (xVal, yVal))

	h.write("    </g>\n")
	h.write("   </g>\n")
	h.write("  </symbol>\n")
	h.write(" </defs>\n")
	h.write("<use xlink:href=\"#chromosome\" x=\"0\" y=\"0\" />\n")
	if count == 2:
		verticalOffset = 400 + spaceBetween - useOuter * 4 * yVal1
		h.write("<use xlink:href=\"#chromosome\" x=\"0\" y=\"0\" transform=\"translate(0 %f) scale(1,-1)\" />\n" % verticalOffset)

def write_svg(dirname, idx, count, frac1, frac2):
	fp.makedirs(dirname)

	with io.open(os.path.join(dirname, "chromosome_%05d.svg" % idx), "w") as h:
		boilerplate_start(h)
		svg_draw(h, count, frac1, frac2)
		boilerplate_end(h)

def main():
	folder = os.path.dirname(__file__)
	svg_dirname = os.path.join(folder, "svg")

	with io.open(os.path.join(folder, "repAnimation.tsv"), "rb") as fh:
		reader = tsv.reader(fh)
		next(reader)
		for line in reader:
			write_svg(svg_dirname, *[float(x) for x in line[:4]])

if __name__ == "__main__":
	main()
