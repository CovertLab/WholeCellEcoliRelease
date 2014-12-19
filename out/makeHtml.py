import os
import re
import argparse

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

def findDirectories(directory):
	allFiles = os.listdir(directory)
	onlyDirs = []
	for f in allFiles:	
		temp = os.path.join(directory,f)
		if os.path.isdir(temp): onlyDirs.append(temp) 
	return onlyDirs

def justName(mystr):
	t = re.split('/',mystr)
	name = t[len(t)-1]
	name = name.replace('.','_')
	return name

def getAlldata(directory):
	data = {'description':'', 'branch': '', 'diff':'', 'hash':''}
	if directory[len(directory)-1] != '/': directory = directory + '/'

	try:
		f = open(directory+'metadata/description')
	except:
		return data
	for line in f: data['description'] = data['description'] + line
	data['description'] = data['description'].strip()
	f.close()

	f = open(os.path.join(directory,'metadata/git_branch'))
	for line in f: data['branch'] = data['branch'] + line
	data['branch'] = data['branch'].strip()
	f.close()

	f = open(os.path.join(directory,'metadata/git_diff'))
	for line in f: data['diff'] = data['diff'] + line
	data['diff'] = data['diff'].strip().replace('\n','<br>').replace('\t',4*'&nbsp;')
	if data['diff'] == '':
		data['diff'] = 'NO DIFF'
	f.close()
	
	f = open(os.path.join(directory,'metadata/git_hash'))
	for line in f: data['hash'] = data['hash'] + line
	data['hash'] = data['hash'].strip()
	f.close()

	return data

def makeHeader(fw, simData):
	fw.write('<html>\n\n')
	fw.write('<head>\n')
	fw.write('	<title>WholeCell E. coli - Output files</title>\n')
	# Format table style
	fw.write('	<style>\n')
	fw.write('	table, th, td {\n')
	fw.write('	    border: 1px solid black;\n')
	fw.write('	    border-collapse: collapse;\n')
	fw.write('	}\n')
	fw.write('	th, td {\n')
	fw.write('	    padding: 15px;\n')
	fw.write('	}\n')
	fw.write('	</style>\n')

	# Format collapsable text style
	fw.write('<style type="text/css">\n')
	fw.write(' .row { vertical-align: top; height:auto !important; }\n')
	fw.write(' .text {display:none; }\n')
	fw.write(' .show {display: none; }\n')
	fw.write(' .hide:target + .show {display: inline; }\n')
	fw.write(' .hide:target {display: none; }\n')
	fw.write(' .hide:target ~ .text {display:inline; }\n')
	fw.write(' @media print { .hide, .show { display: none; } }\n')
	fw.write(' </style>\n')

	fw.write('</head>\n\n')

	fw.write('<script>\n\n')
	
	for i in simData:
		getDescriptionData = getAlldata(i) 
		cont = 'description:'+ getDescriptionData['description']+'<br>branch: '+getDescriptionData['branch']+'<br>'
		for j in simData[i]:
			namei = justName(i)
			fw.write('function getUrl_'+namei+ '_'+j+'()\n')
			fw.write('{\n')
	 		#all img file names
	 		fw.write('	var fileNames = ["'+simData[i][j]['files'][0]+'"')
	 		for k in range(1,len(simData[i][j]['files'])):
	 			fw.write(',"'+simData[i][j]['files'][k]+'"')
	 		fw.write('];\n')
	 		fw.write('	var directory = "'+simData[i][j]['dir']+'";\n')
	 		fw.write('	var contents= "'+cont+'";\n')

	 		fw.write('	directory = encodeURIComponent(directory);\n')
	 		fw.write('	fileNames = encodeURIComponent(fileNames);\n')
	 		fw.write('	contents = encodeURIComponent(contents);\n')

			fw.write('	return \"images.html?var1=\" + fileNames + \"&dir=\" + directory + \"&content=\" + contents;\n')
			fw.write('}\n\n')
    
	fw.write('</script>\n\n')

def makeBody(fw, simData):
	fw.write('<body>\n')

	fw.write('<table style="width:100%">\n')
	fw.write('  <tr>\n')
	fw.write('    <td>Description</td>\n')
	fw.write('    <td>Branch</td> \n')
	fw.write('    <td>Difference</td>\n')
	fw.write('    <td>Hash</td>\n')
	fw.write('    <td>Simulations</td>\n')
	fw.write('  </tr>\n')

	for idx,i in enumerate(simData):
		if len(simData[i]) == 0: continue
		getDescriptionData = getAlldata(i) 
		fw.write('  <tr>\n')
		fw.write('    <td>'+getDescriptionData['description']+'</td>\n')
		fw.write('    <td>'+getDescriptionData['branch']+'</td> \n')

		# Create expand/collapse function for diff
		fw.write('		<td>\n')
		fw.write('		<div class="row">\n')
		fw.write('		 <a href="#hide' + str(idx) + '" class="hide" id="hide' + str(idx) + '">Expand</a>\n')
		fw.write('		 <a href="#show' + str(idx) + '" class="show" id="show' + str(idx) + '">Collapse</a>\n')
		fw.write('		 <div class="text">\n')
		fw.write('		 <font size="1" face="Courier"><br>')
		fw.write('			'+getDescriptionData['diff']+'\n')
		fw.write('		 </font>')
		fw.write('		</div>\n</div>\n')
		fw.write('		</td>\n')

		fw.write('    <td>'+getDescriptionData['hash']+'</td>\n')
		fw.write('    <td>')
		for j in simData[i]:
			namei = justName(i)
			url = 'getUrl_'+namei+ '_'+j+'()'
			fw.write('        <a href=\"javascript:document.location.href='+url+';\" target=\"_blank\">'+namei+'_'+j+'</a><br>\n')
		fw.write('    </td>\n')
		fw.write('  </tr>\n')

	fw.write('</table>\n')
	fw.write('</body>\n')
	fw.write('</html>\n')


def main(out_directory):
		
	outDirectory = out_directory
	
	allSimulations = findDirectories(outDirectory)

	allSimulationsData = {}	
	for i in allSimulations:
		dirs = findDirectories(i)
		if dirs == []: continue
		allSimulationsData[i] = {}
		for j in dirs:
			files = findFiles(os.path.join(j,'plotOut'), '.svg')
			if files == []: continue
			myDir = j
			if '/' != j[len(j)-1]: myDir = j + '/'
			myDir = myDir.replace('home/users','Volumes')
			namej = justName(j)
			allSimulationsData[i][namej] = {'files': files, 'dir': myDir}


	#make htmlfile 
	if '/' != outDirectory[len(outDirectory)-1]: outDirectory = outDirectory + '/'
	htmlFile = outDirectory+'index.html'
	fw = open(htmlFile, 'w')
	makeHeader(fw, allSimulationsData)
	makeBody(fw, allSimulationsData)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("out_directory", type = str)

	args = parser.parse_args().__dict__

	main(args["out_directory"])
