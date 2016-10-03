import sys
import os
import glob
import re

def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

nrOfDirs = len(sys.argv) - 1

for i in range(1,(nrOfDirs + 1)):
	path = sys.argv[i]
	print os.path.basename(os.path.normpath(path))
	outputdir = []
	
	for files in glob.glob(path + '*.txt'):
		name = os.path.basename(os.path.normpath(files))
		name = name.split(".")[0]
		IN = open(files, 'r')
		for j in range(0,48):
			IN.next()
		averageAllowedDev = IN.next().rstrip().split(":")[1]
		averageAllowedDev = averageAllowedDev[1:6]
		
		binTest = []
		for line in IN:
			#print line[0:26]
			if str(line[0:19]) == "Windowed, bin test:":
				break
		IN.next()
		binTestLine = IN.next()
		if str(binTestLine)[0:6] != "Found:":
			#print "!!!"
			#print str(binTestLine)[0:6]
			binTestLine = binTestLine.rstrip().lstrip().split("\t")
			binTest.append(binTestLine)
			testing = True
			while testing == True:
				binTestLine = IN.next()
				binTestLine = binTestLine.rstrip().lstrip().split("\t")
				#print "tweede?"
				#print binTestLine
				if binTestLine[0] == "+" or binTestLine[0] == "-":
					#print "ja!"
					binTest.append(binTestLine)
				else:
					#print "nee"
					testing = False
			
		for line in IN:
			if str(line[0:26]) == "Windowed, aneuploidy test:":
				break
		aneu = IN.next().lstrip().rstrip()
		
		if str(aneu) != "Nothing found":
			outputdir.append(name + "\t" + aneu + "\t" + averageAllowedDev)
		elif binTest:
			teststring = []
			for tests in binTest:
				#print tests
				teststring.append(tests[4].split(":")[0] + ", " + tests[3] + ", (" + tests[0] + ")")
			outputdir.append(name + "\t" + " & ".join(teststring) + "\t" + averageAllowedDev)
			
		else:
			outputdir.append(name + "\t0\t" + averageAllowedDev)
	
	outputdir.sort(key=natural_key)		
	print "\n".join(outputdir)
		
			
