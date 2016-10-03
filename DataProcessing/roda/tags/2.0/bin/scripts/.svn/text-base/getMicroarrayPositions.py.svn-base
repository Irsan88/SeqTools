##############################################################################
# Usage: python getMicroarrayPositions.py [array cyhd.txt] [output.bed]      #
# [positions filter.bed]                                                     #
#                                                                            #
# This script exports the microarray positions as .bed file. This file can   #
# be used for making the samtools mpileup file which is required for the     #
# SNP concordance check.                                                     #
#                                                                            #
##############################################################################

import sys

arrayFile = sys.argv[1]
outputFile = sys.argv[2]
bedFile = sys.argv[3]

arraySNPs = []



def loadBed(inputFile):
	listBed = open(inputFile, 'r')
	arrayBed = {}
	listBed.next()			# Skip header line
	for lines in listBed:
		linewords = lines.rstrip().split('\t')
		if linewords[0][3:] not in arrayBed:
			arrayBed[linewords[0][3:]] = []
		arrayBed[linewords[0][3:]].append([int(linewords[1]), int(linewords[2])])
	return arrayBed

def liesInBed(chrom, position, arrayBed):
	#print "start"
	#print chrom
	if chrom in arrayBed:
		for regions in arrayBed[chrom]:
			if position >= regions[0] and position <= regions[1]:
				#print "lies in bed!"
				#print regions
				return True
			elif position > regions[1]:
				#print "bij continue"
				continue
			elif position < regions[0]:
				#print "bij laatste"
				#print regions
				#print position
				return False
	else:
		#print "hier helemaal"
		return False

def main():
	arrayBed = loadBed(bedFile)

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
				if linewords[7] not in "M" and liesInBed(linewords[7], int(linewords[8]), arrayBed): #and float(linewords[2]) > 0.01:
					arraySNPs.append([linewords[7], linewords[8]])
	inFile.close()

	outFile = open(outputFile, 'w')
	for positions in arraySNPs:
		outFile.write("chr" + '\t'.join(positions) + "\n")

	outFile.close()

if __name__ == "__main__":
    main()

