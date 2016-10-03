##############################################################################
# Usage: python snpConcordance.py [array cyhd.txt] [NGS snps.filt.vcf]       #
# [depthfile.txt]                                                            #
#                                                                            #
# This script compares the SNPs called by both SNP array and NGS. The counts # 
# are returned, as well as the corresponding SNPs (unless printing SNPs is   #
# disabled with 'Y' at argument 3).                                          #
#                                                                            #
# Testing resulted in following:                                             #
# % of array/NGS calls matching the position of the other: 
# % of position matches that have same genotype for array/NGS:
# % of position matches that do not have same genotype for array/NGS:
#
##############################################################################

import sys

usage = "Wrong input!\nUsage:\tpython snpConcordance.py [arrayfile.cyhd.txt] [NGS variants.vcf] [targeted depth file] [positions array filtered by bed] [output.concor]"

if len(sys.argv) != 6:
	print usage
	sys.exit()

filterConfidence = 0.01
arrayFile = sys.argv[1]
ngsFile = sys.argv[2]
depthFile = sys.argv[3]
positionsFilteredBed = sys.argv[4]
outputFile = sys.argv[5]
printReads = False

def filterOnBed(inputdict, bedPositions):
	outputdict = {}
	for chrom in inputdict:
		for snp in inputdict[chrom]:
			if chrom in bedPositions:
				if int(snp[0]) in bedPositions[chrom]:
					if chrom not in outputdict:
						outputdict[chrom] = []
					outputdict[chrom].append(snp)
	return outputdict
	
def checkForReadEnds(string):
	index = string.find('^')
	while index >= 0:
		string = string[:index] + string[index + 2:]
		index = string.find('^')
	return string
	
def getGenotype(mpileup):
	threshold = 0.2
	ref = mpileup[1].upper()
	cov = int(mpileup[2])
	calls = ''.join(c for c in mpileup[3] if c not in '$')
	calls = checkForReadEnds(calls)
	numberOfReads = len(calls)
	refFreq = sum([x in ',.' for x in calls]) / float(numberOfReads)
	#print refFreq
	if refFreq >= 1 - threshold:
		return "0/0", ref, ref
	freqs = {}
	for char in calls:
		if char not in ',.':
			char = char.upper()
			if char not in freqs:
				freqs[char] = 0
			freqs[char] += 1
	for char in freqs:
		freqs[char] = freqs[char] / float(numberOfReads)
	freqs['ref'] = refFreq
	#print freqs		
	
	alt = []
	for char in freqs:
		if freqs[char] >= 1 - threshold:
			return "1/1", char, char
		elif refFreq >= threshold and char >= threshold:
			return "0/1", ref, char
		elif char >= threshold:
			alt.append(char)
	if len(alt) == 2:
		return "1/2", alt[0], alt[1]
	else:
		OUT.write("Vreemd!\t" + str(calls) + "\n" + str(freqs) + "\n")
		
	return None, None, None
	

bedPositions = {}
inFile = open(positionsFilteredBed, 'r')
for line in inFile:
	line = line.rstrip().split("\t")
	line[0] = line[0][3:]
	if line[0] not in bedPositions:
		bedPositions[line[0]] = []
	bedPositions[line[0]].append(int(line[1]))

arraySNPs = {}
ngsSNPs = {}
depthSNPs = {}

countje = 0
qualityFiltered = {}
inFile = open(arrayFile, 'r')
lines = inFile.readlines()
first = True
for line in lines:
	if line[0] not in '#':
		linewords = line.rstrip().split('\t')
		if first:
			arHead = [linewords[0], linewords[8], linewords[6], linewords[1], linewords[5], linewords[2], linewords[3], linewords[4], linewords[0]]
			first = False
		else:
			if linewords[7] not in "M" and float(linewords[2]) <= filterConfidence:
				countje += 1
				#if liesInBed(linewords[7], int(linewords[8])):
				if linewords[7] not in arraySNPs:
					arraySNPs[linewords[7]] = []
				arraySNPs[linewords[7]].append([linewords[8], linewords[6], linewords[1], linewords[5], linewords[2], linewords[3], linewords[4], linewords[0]])
			elif linewords[7] not in "M" and float(linewords[2]) > filterConfidence:
				if linewords[7] not in qualityFiltered:
					qualityFiltered[linewords[7]] = []
				qualityFiltered[linewords[7]].append(linewords[8])
inFile.close()
#print "Before filtering on bed:\t" + str(countje)
arraySNPsFilt = filterOnBed(arraySNPs, bedPositions)
#print "Array read"		

'''
count2 = 0
count3 = 0
for chrom in arraySNPsFilt:
	for snp in arraySNPsFilt[chrom]:
		count2 += 1
		if float(snp[4]) >= filterConfidence:
			count3 += 1
print "After filtering on bed:\t" + str(count2)
print "After filtering on quality:\t" + str(count3)
'''			


inFile = open(ngsFile, 'r')
lines = inFile.readlines()
for line in lines:
	if line[0:2] not in '##':
		linewords = line.rstrip().split('\t')
		if line[0] == "#":
			ngsHead = linewords
		else:
			linewords[0] = linewords[0][3:]
			if linewords[0] not in "M" and linewords[6] == "PASS":
				if linewords[0] in qualityFiltered:
					if linewords[1] in qualityFiltered[linewords[0]]:
						continue
				if linewords[0] not in ngsSNPs:
					ngsSNPs[linewords[0]] = []
				ngsSNPs[linewords[0]].append(linewords[1:])
inFile.close()
ngsSNPsFilt = filterOnBed(ngsSNPs, bedPositions)
#print "NGS read"

inFile = open(depthFile, 'r')
lines = inFile.readlines()
for line in lines:
	linewords = line.rstrip().split('\t')
	linewords[0] = linewords[0][3:]
	if linewords[0] in qualityFiltered:
		if linewords[1] in qualityFiltered[linewords[0]]:
			continue
	if linewords[0] not in depthSNPs:
		depthSNPs[linewords[0]] = []
	depthSNPs[linewords[0]].append(linewords[1:])
inFile.close()
#print "Depthfile read"

countNgsTotal = 0
countArrayTotal = 0
countDepthTotal = 0
for chrom in ngsSNPsFilt:
	countNgsTotal = countNgsTotal + len(ngsSNPsFilt[chrom])
for chrom in arraySNPsFilt:
	countArrayTotal = countArrayTotal + len(arraySNPsFilt[chrom])
for chrom in depthSNPs:
	countDepthTotal = countDepthTotal + len(depthSNPs[chrom])

OUT = open(outputFile, 'w')

OUT.write("\n##SNPs total NGS:\t" + str(countNgsTotal) + "\n")
OUT.write("##SNPs total Array:\t" + str(countArrayTotal) + "\n")
OUT.write("##SNPs total Depth:\t" + str(countDepthTotal) + "\n")

####################################
# Step 1 - Compare NGS variantlist #
####################################

arrayPos = {}
ngsPos = {}
depthPos = {}

for chrom in arraySNPsFilt:	
	arrayPos[chrom] = [snp[0] for snp in arraySNPsFilt[chrom]]
for chrom in ngsSNPsFilt:
	ngsPos[chrom] = [snp[0] for snp in ngsSNPsFilt[chrom]]	
for chrom in depthSNPs:
	depthPos[chrom] = [snp[0] for snp in depthSNPs[chrom]]

#Search for array SNPs in NGS
arrayRef = {}
for chrom in arraySNPsFilt:
	for snp in arraySNPsFilt[chrom]: 
		if chrom in ngsPos:
			if snp[0] in ngsPos[chrom]:
				if chrom not in arrayRef:
					arrayRef[chrom] = []
				arrayRef[chrom].append(snp)
				
	
#Search for NGS SNPs from VCF in array
ngsRef = {}
for chrom in ngsSNPsFilt:	#ngsResult:
	for snp in ngsSNPsFilt[chrom]:	
		if chrom in arrayPos:
			if snp[0] in arrayPos[chrom]:
				if chrom not in ngsRef:
					ngsRef[chrom] = []
				ngsRef[chrom].append(snp)
	
#Delete SNPs found in both from arraySNPs
for chrom in arrayRef:
	for snp in arrayRef[chrom]:
		arraySNPsFilt[chrom].remove(snp)	

#Delete SNPs found in both from ngsSNPs
for chrom in ngsRef:
	for snp in ngsRef[chrom]:
		ngsSNPsFilt[chrom].remove(snp)	


#################################
# Step 2 - Compare depth file   #
#################################
#Here we remove the NGS variants from 
ngsRefPos = {}
for chrom in ngsRef:
	ngsRefPos[chrom] = [snp[0] for snp in ngsRef[chrom]]

deleteFromDepth = {}
for chrom in depthSNPs:
	for snp in depthSNPs[chrom]:
		if chrom in ngsRefPos:
			if snp[0] in ngsRefPos[chrom]:
				if chrom not in deleteFromDepth:
					deleteFromDepth[chrom] = []
				deleteFromDepth[chrom].append(snp)

for chrom in deleteFromDepth:
	for snp in deleteFromDepth[chrom]:
		depthSNPs[chrom].remove(snp)

depthPos = {}
for chrom in depthSNPs:
	depthPos[chrom] = [snp[0] for snp in depthSNPs[chrom]]

arrayPos = {}
for chrom in arraySNPsFilt:
	arrayPos[chrom] = [snp[0] for snp in arraySNPsFilt[chrom]]	

#Search for remaining positions of the array in the depth file
depthRef = {}
for chrom in depthSNPs:
	for snp in depthSNPs[chrom]:
		if chrom in arrayPos:
			if snp[0] in arrayPos[chrom]:
				if chrom not in depthRef:
					depthRef[chrom] = []
				depthRef[chrom].append(snp)
				
#Search for remaining positions of the depthfile in the array
arrayRef2 = {}
for chrom in arraySNPsFilt:
	for snp in arraySNPsFilt[chrom]: 
		if chrom in depthPos:
			if snp[0] in depthPos[chrom]:
				if chrom not in arrayRef2:
					arrayRef2[chrom] = []
				arrayRef2[chrom].append(snp)


#Delete SNPs found in both from arraySNPsFilt
for chrom in arrayRef2:
	for snp in arrayRef2[chrom]:
			arraySNPsFilt[chrom].remove(snp)	

#Delete SNPs found in both from depthSNPs
for chrom in depthRef:
	for snp in depthRef[chrom]:
			depthSNPs[chrom].remove(snp)	


####################################################
# Step 1 - Array and NGS variants genotype compare #
####################################################

ngsNoMatch = {}
ngsMatch = {}
for chrom in ngsRef:
	#print "Chromosome " + str(chrom)
	for snp in ngsRef[chrom]:
		arrayCheck = ""
		genoType = snp[-1].split(":")[0]
		ref = snp[2]
		alt = snp[3]
		if genoType in "1/1":
			ref = alt
		elif genoType in "0/0":
			alt = ref
		elif genoType in "0/1":
			pass
		elif genoType in "1/2":
			alt = alt.split(",")
			if len(alt) > 2:
				OUT.write("More than two alleles: " + str(alt) + "\n")
			ref = alt[0]
			alt = alt[1]
		else:
			OUT.write(str(genoType) + "\n")
		#print "SNP done: " + str(snp) + "\n" + str(alt+ref) + "\n" + str(genoType)
		for index, value in enumerate(arrayRef[chrom]):
			if value[0] == snp[0]:
				arrayCheck = arrayRef[chrom][index]
				if arrayCheck[3] in str(ref+alt) or arrayCheck[3] in str(alt+ref):
					if chrom not in ngsMatch:
						ngsMatch[chrom] = []
					ngsMatch[chrom].append([snp, arrayCheck])
					#ngsMatch[chrom].append(arrayCheck)
				else:
					if chrom not in ngsNoMatch:
						ngsNoMatch[chrom] = []
					ngsNoMatch[chrom].append([snp, arrayCheck])
					#ngsNoMatch[chrom].append(arrayCheck)
				break

###################################
#See if the observed genotype of depthRef and arrayRef also matches at the positions
######################################
depthNoMatch = {}
depthMatch = {}
for chrom in depthRef:
	#print "Chromosome " + str(chrom)
	for snp in depthRef[chrom]:
		arrayCheck = ""
		(genoType, ref, alt) = getGenotype(snp)
		#print "SNP done: " + str(snp) + "\n" + str(alt) + str(ref) + "\n" + str(genoType)
		for index, value in enumerate(arrayRef2[chrom]):
			if value[0] == snp[0]:
				arrayCheck = arrayRef2[chrom][index]
				if arrayCheck[3] in str(ref+alt) or arrayCheck[3] in str(alt+ref):
					if chrom not in depthMatch:
						depthMatch[chrom] = []
					depthMatch[chrom].append([snp, arrayCheck])
					#depthMatch[chrom].append(arrayCheck)
				else:
					if chrom not in depthNoMatch:
						depthNoMatch[chrom] = []
					depthNoMatch[chrom].append([snp, arrayCheck])
					#depthNoMatch[chrom].append(arrayCheck)
				break
		
countArraySNPs = 0
for chrom in arraySNPsFilt:
	countArraySNPs = countArraySNPs + len(arraySNPsFilt[chrom])
OUT.write('\n##SNPs only in array:\t' + str(countArraySNPs) + "\n")
if printReads:
	for chrom in sorted(arraySNPsFilt):
		for snp in arraySNPsFilt[chrom]:
			OUT.write(str(chrom) + "\t" + '\t'.join(snp) + "\n")

countNgsSNPs = 0
for chrom in ngsSNPsFilt:
	countNgsSNPs = countNgsSNPs + len(ngsSNPsFilt[chrom])
OUT.write('##SNPs only in ngs variantlist:\t' + str(countNgsSNPs) + "\n")
if printReads:
	for chrom in sorted(ngsSNPsFilt):
		for snp in ngsSNPsFilt[chrom]:
			OUT.write(str(chrom) + "\t" + '\t'.join(snp) + "\n")
			
countDepthSNPs = 0
for chrom in depthSNPs:
	countDepthSNPs = countDepthSNPs + len(depthSNPs[chrom])
OUT.write('##SNPs only in ngs depthfile:\t' + str(countDepthSNPs) + "\n")
if printReads:
	for chrom in sorted(depthSNPs):
		for snp in depthSNPs[chrom]:
			OUT.write(str(chrom) + "\t" + '\t'.join(snp) + "\n")

countNGS = 0
countArray = 0
for chrom in arrayRef:
	countArray = countArray + len(arrayRef[chrom])
for chrom in ngsRef:
	countNGS = countNGS + len(ngsRef[chrom])
OUT.write('\n##SNPs in both variantlist and array:\t' + str(len(arrayRef.keys())) + " chromosomes,\t" + str(countArray) + " SNPs" + "\n")

countDepth = 0
countArray2 = 0
for chrom in arrayRef2:
	countArray2 = countArray2 + len(arrayRef2[chrom])
for chrom in depthRef:
	countDepth = countDepth + len(depthRef[chrom])
OUT.write('##SNPs in both depthfile and array:\t' + str(len(arrayRef2.keys())) + " chromosomes,\t" + str(countArray2) + " SNPs" + "\n")


countMatch = 0
countNoMatch = 0
for chrom in ngsMatch:
	countMatch = countMatch + len(ngsMatch[chrom])
for chrom in ngsNoMatch:
	countNoMatch = countNoMatch + len(ngsNoMatch[chrom])

OUT.write('\n##Variants - Array No match:\t' + str(countNoMatch) + "\n")
if printReads:
	for chrom in sorted(ngsNoMatch):
		for snp in ngsNoMatch[chrom]:
			OUT.write("\nNGS SNP:\t" + str(chrom) + "\t" + str(snp[0]) + "\n")
			OUT.write("Array SNP:\t" + str(chrom) + "\t" + str(snp[1]) + "\n")

OUT.write('##Variants - Array Match:\t' + str(countMatch) + "\n")
if printReads:
	for chrom in sorted(ngsMatch):
		for snp in ngsMatch[chrom]:
			OUT.write("\nNGS SNP:\t" + str(chrom) + "\t" + str(snp[0]) + "\n")
			OUT.write("Array SNP:\t" + str(chrom) + "\t" + str(snp[1]) + "\n")


countDepthMatch = 0
countDepthNoMatch = 0
for chrom in depthMatch:
	countDepthMatch = countDepthMatch + len(depthMatch[chrom])
for chrom in depthNoMatch:
	countDepthNoMatch = countDepthNoMatch + len(depthNoMatch[chrom])

OUT.write('##Depthfile - Array No match:\t' + str(countDepthNoMatch) + "\n")
if printReads:
	for chrom in sorted(depthNoMatch):
		for snp in depthNoMatch[chrom]:
			pOUT.write("\nNGS SNP:\t" + str(chrom) + "\t" + str(snp[0]) + "\n")
			OUT.write("Array SNP:\t" + str(chrom) + "\t" + str(snp[1]) + "\n")

OUT.write('##Depthfile - Array Match:\t' + str(countDepthMatch) + "\n")
if printReads:
	for chrom in sorted(depthMatch):
		for snp in depthMatch[chrom]:
			OUT.write("\nNGS SNP:\t" + str(chrom) + "\t" + str(snp[0]) + "\n")
			OUT.write("Array SNP:\t" + str(chrom) + "\t" + str(snp[1]) + "\n")

OUT.write('\n##Total SNPs that do not match:\t' + str(countNoMatch + countDepthNoMatch) + "\n")
OUT.write('##Total SNPs that match:\t' + str(countMatch + countDepthMatch) + "\n")

totalCalls = countNoMatch + countDepthNoMatch + countMatch + countDepthMatch + countDepthSNPs + countNgsSNPs + countArraySNPs
OUT.write('\n\nSNP CONCORDANCE:\t' + str(round((countMatch + countDepthMatch) / float(totalCalls) * 100, 2)) + str(" %") + "\n")

OUT.close()
