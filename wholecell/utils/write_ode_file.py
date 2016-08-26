import numpy as np

def writeOdeFile(fileName, derivatives, derivativesJacobian):
	h = open(fileName, "w")
	h.write("import numpy as np\n")
	h.write("\n")
	h.write("def derivatives(y, t):\n")
	h.write("\treturn np.array(" + str(derivatives)[7:-1] + ").reshape(-1)\n")
	h.write("\n")
	h.write("def derivativesJacobian(y, t):\n")
	h.write("\treturn np.array(" + str(derivativesJacobian)[7:-1] + ")\n")
	h.close()