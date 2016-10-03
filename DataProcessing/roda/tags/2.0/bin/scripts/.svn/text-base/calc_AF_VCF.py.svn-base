import sys

inputFile = sys.argv[1]
outputFile = sys.argv[2]

IN = open(inputFile, 'r')
OUT = open(outputFile, 'w')

for line in IN:
	if line[0] == '#':
		OUT.write(line)
	else:
		line = line.rstrip().split('\t')
		if line[8] == 'GT:AD:DP:GQ:PL':
			(GT, AD, DP, GQ, PL) = line[9].split(':')
		elif line[8] == 'GT:AD:GQ:PL':
			(GT, AD, GQ, PL) = line[9].split(':')
		else:
			print '\t'.join(line) + "\nFormat of the column differs"
			continue
		freqs = []
		AD = AD.split(',')
		totalReads = 0
		for allele in AD:
			totalReads = totalReads + int(allele)
			freqs.append(int(allele))
		newColumn = [""]
		line[7] = line[7].split(';')
		for entry in line[7]:
			if entry[:2] != "AF":
				newColumn[-1] = newColumn[-1] + str(entry) + ";"
			else:
				newColumn.append("")
		if totalReads != 0:
			freqs = [('%f' % (round(x/float(totalReads), 3))).rstrip('0').rstrip('.') for x in freqs]
			OUT.write('\t'.join(line[0:7]) + '\t' + newColumn[0] + 'AF=' + ','.join(freqs[1:]) + ';' + newColumn[1] + '\t' + '\t'.join(line[8:10]) + '\n')
		else:
			noAlleles = len(freqs) - 1
			newFreqs = [0] * noAlleles
			newFreqs = [str(x) for x in newFreqs]
			OUT.write('\t'.join(line[0:7]) + '\t' + newColumn[0] + 'AF=' + ','.join(newFreqs) + ';' + newColumn[1] + '\t' + '\t'.join(line[8:10]) + '\n')
	

IN.close()
OUT.close()		
