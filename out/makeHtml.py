import os
import re
from collections import OrderedDict
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
		if f == 'kb' or f == 'metadata': continue
		temp = os.path.join(directory,f)
		if os.path.isdir(temp): onlyDirs.append(temp) 
	onlyDirs.sort()
	onlyDirs.reverse()
	return onlyDirs

def justName(mystr):
	t = re.split('/',mystr)
	name = t[len(t)-1]
	name = name.replace('.','_')
	return name

def getAlldata(directory, flag):
	data = {'description':'', 'branch': '', 'diff':'', 'hash':'', 'short_name': ''}
	if directory[len(directory)-1] != '/': directory = directory + '/'

	try:
		f = open(directory+'metadata/description')
	except:
		return data
	for line in f: data['description'] = data['description'] + line
	data['description'] = data['description'].strip()
	f.close()

	if flag:
		try:
			f = open(directory+'metadata/short_name')
		except:
			return data
		for line in f: data['short_name'] = data['short_name'] + line
		data['short_name'] = data['short_name'].strip()
		f.close()

		return data

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
		for k in simData[i]:
			getDescriptionData = getAlldata(i+'/'+k,1) 
			cont = 'description: '+ getDescriptionData['description']+'<br>short_name: '+getDescriptionData['short_name']+'<br>'

			for j in simData[i][k]:
				namei = justName(i)
				namek = justName(k)
				fw.write('function getUrl_'+namei+ '_'+namek+'_'+j+'()\n')
				fw.write('{\n')
		 		#all img file names
		 		fw.write('	var fileNames = ["'+simData[i][k][j]['files'][0]+'"')
		 		for l in range(1,len(simData[i][k][j]['files'])):
		 			fw.write(',"'+simData[i][k][j]['files'][l]+'"')
		 		fw.write('];\n')
		 		fw.write('	var directory = "'+simData[i][k][j]['dir']+'";\n')
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
		getDescriptionData = getAlldata(i,0) 
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

		fw.write('        <table style="width:100%">\n')
		fw.write('        <tr>\n')
		fw.write('            <td>Description</td>\n')
		fw.write('            <td>short_name</td> \n')
		fw.write('            <td>Links</td>\n')
		fw.write('        </tr>\n')
	
		for k in simData[i]:
			desc = getAlldata(i+'/'+k,1)
			namek = justName(k)
			fw.write('        <tr>\n')
			fw.write('            <td>'+desc['description']+'</td>\n')
			fw.write('            <td>'+desc['short_name']+'</td> \n')
	
			fw.write('    		  <td>')
			for j in simData[i][k]:
				namei = justName(i)
				url = 'getUrl_'+namei+ '_'+namek+ '_'+j+'()'
				fw.write('        <a href=\"javascript:document.location.href='+url+';\" target=\"_blank\">'+namei+'_'+namek+'_'+j+'</a><br>\n')
			fw.write('            </td>\n')
			fw.write('        <tr>\n')
		fw.write('        </table>\n')

		fw.write('    </td>\n')
		fw.write('  </tr>\n')

	fw.write('</table>\n')
	fw.write('</body>\n')
	fw.write('</html>\n')


def main(out_directory):
		
	outDirectory = out_directory
	
	allSimulations = findDirectories(outDirectory)

	allSimulationsData = OrderedDict({})	
	for i in allSimulations:
		dirs = findDirectories(i)
		if dirs == []: continue
		allSimulationsData[i] = {}
		for j in dirs:
			subDirs = findDirectories(j)
			if subDirs == []: continue
			allData = OrderedDict({})
			for k in subDirs:
				files = findFiles(os.path.join(k,'plotOut'), '.svg')
				if files == []: continue
				
				myDir = k
				if '/' != k[len(k)-1]: myDir = k + '/'
				myDir = myDir.replace('home/users','Volumes')
				namek = justName(k)
				allData[namek] = {'files': files, 'dir': myDir}
			namej = justName(j)
			allSimulationsData[i][namej] = allData

	#make htmlfile 
	if '/' != outDirectory[len(outDirectory)-1]: outDirectory = outDirectory + '/'
	htmlFile = outDirectory+'index.html'
	import ipdb; ipdb.set_trace()
	fw = open(htmlFile, 'w')
	makeHeader(fw, allSimulationsData)
	makeBody(fw, allSimulationsData)
	fw.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("out_directory", type = str)

	args = parser.parse_args().__dict__

	main(args["out_directory"])
