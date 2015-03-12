from pylint import epylint as lint

def writeOnFile(myBuffer, filename):
	f = open(filename,'w')
	for line in myBuffer: 
		f.write(line)
	f.close()

def checkPythonModule(moduleName, options, outputErrorFile, outputWarningFile):
	(pylint_stdout, pylint_stderr) = lint.py_run(moduleName+' '+options, return_std=True)
	print 'Done py_run()'
	writeOnFile(pylint_stdout, outputWarningFile)
	writeOnFile(pylint_stderr, outputErrorFile)


#call functions
outputErrorFile = '/home/users/sajia/wcEcoli/pylintError.txt'
outputWarningFile = '/home/users/sajia/wcEcoli/pylintWarning.txt'
moduleName = '/home/users/sajia/wcEcoli/wholecell/'
options = '-d W0312,C0111,C0330'

checkPythonModule(moduleName, options, outputErrorFile, outputWarningFile)
 
#from command line: pylint -d W0312,C0111,C0330 /home/users/sajia/wcEcoli/wholecell/ > outputfile.txt