import numpy as np

def writeMetabolicConstraintsFile(fileName, constraints):
	h = open(fileName, "w")
	h.write("import numpy as np\n")
	h.write("\n")
	h.write("def constraints(enzymes, kineticsSubstrates):\n")
	# If you want to make the file readable (with lots of line breaks) change the "str" to "repr" in the line below
	h.write("\treturn np.array(" + str(constraints)[7:-1] + ").reshape(-1)\n")
	h.write("\n")
	h.close()
