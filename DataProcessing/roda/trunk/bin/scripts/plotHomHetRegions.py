import sys
import matplotlib.pyplot as plt

vcfIn = open(sys.argv[1], 'r')

positions = {}
allChroms = []
lengthContigs = {}

for line in vcfIn:
	if line[0] == '#':
		if line[0:8] == '##contig':
			line = line.rstrip().split(',')
			chrom = line[0].split('=')[-1]
			length = line[1].split('=')[-1]
			lengthContigs[chrom] = length
		continue
	line = line.rstrip().split('\t')
	if line[0] not in positions:
		positions[line[0]] = []
		allChroms.append(line[0])
	alleleFreq = line[7].split(';')[1]
	alleleFreq = alleleFreq.split('=')
	if alleleFreq[0] != 'AF':
		print "raar!" + str(alleleFreq)
	alleleFreq = alleleFreq[1]
	positions[line[0]].append([line[1], alleleFreq])

#subplots = [0, 611, 612, 613, 614, 615, 616]
#subplots = [0, 511, 512, 513, 514, 515]
subplots = [0, 411, 412, 413, 414]
loop = 0
nrFig = 0
allChroms.sort(key=lambda x: x[0])
fig = plt.figure()

for chrom in allChroms:
	loop = loop + 1
	if loop == 5:
		nrFig = nrFig + 1
		plt.savefig(sys.argv[2] + "_" + str(nrFig) + '.png')
		fig = plt.figure()
		loop = 1
	print chrom
	posx = []
	posy = []
	for pos in positions[chrom]:
		freqs = pos[1].split(',')
		lengthfreqs = len(freqs)
		for i in range(0,lengthfreqs):
			posx.append(int(pos[0]))
			posy.append(float(freqs[i]))	
		
	
	ax = fig.add_subplot(subplots[loop])
	ax.scatter(posx, posy, s=2, facecolor='0.5', lw = 0)
	ax.set_ylabel('AF ' + str(chrom))
	#ax.set_xlabel('Position on ' + str(chrom))
	ax.set_ylim(-0.1,1.1)
	ax.set_xlim(0, int(lengthContigs[chrom]))
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
		item.set_fontsize(10)
	

nrFig = nrFig + 1
plt.savefig(sys.argv[2] + "_" + str(nrFig) + '.png')
