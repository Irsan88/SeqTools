#!/usr/bin/env python

# Script based on qc.R R-script by sdentro
# Daphne van Beek

import matplotlib
import sys
from pylab import *
from reportlab.lib import colors as ncolors
from reportlab.lib.pagesizes import letter, inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Image, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
import os
import argparse
import glob
	
parser = argparse.ArgumentParser(description='Create .bed file for exon kit')
parser.add_argument('pathSample', metavar='[path to sample without extension]', type=str, help='path to sample without extension')
parser.add_argument('-generalbed', metavar='[general.bed]', type=str, help='path to general.bed file')
parser.add_argument('-exon', metavar='[exon start and stop.bed]', type=str, help='path to exon start and stop.bed file')
parser.add_argument('-panel', metavar='[path to genepanel without extension]', type=str, help='path to genepanel file (.._genePanels) without extension')
parser.add_argument('-hist', action='store_true', help='Y if hist file exists.')

f = parser.parse_args()	
	
#############
# Set paths	#
#############
generalBed = None 
path_to_snps_metric = None
path_to_config = glob.glob(f.pathSample + ".conf*")
if f.generalbed:
	generalBed = f.generalbed
bedIn = None
if f.exon:
	bedIn = f.exon
	path_to_snps_metric = f.pathSample + ".stat.snps.metrics"		#".snps.metrics"
path_to_flagstat = f.pathSample + ".stat.mark.bam.flagstat"#".stat.recal.bam.flagstat"		#".recal.bam.flagstat"
if bedIn:
	depthFile = f.pathSample + ".stat.recal.bam.depth"			#".recal.bam.depth"	
if bedIn or f.hist:
	path_to_target_flagstat = f.pathSample + ".stat.target.recal.bam.flagstat"		#".target.recal.bam.flagstat"
if f.hist:
	path_to_histFile = f.pathSample + ".stat.recal.bam.hist"
outputPDF = SimpleDocTemplate(f.pathSample + '.pdf', pagesize=letter)
sampleName = str(f.pathSample.split('/')[-1])


#########################
# Set (style) defaults	#
#########################
minimum_depth = 30
aimed_depth = 50
depth_threshold = 30
styles = getSampleStyleSheet()
styleHeading = styles['Heading2']
styleHeading.alignment = 1
styleNormal = styles['Heading4']
styleNormal.alignment = 1
elements = []

Ncols = 3
figheight = 4.4   # inches
figwidth = figheight * Ncols    # inches
figure(1, figsize=(figwidth, figheight))
rcParams['font.size'] = 10.0
rcParams['axes.titlesize'] = 16.0
rcParams['xtick.labelsize'] = 12.0
rcParams['legend.fontsize'] = 12.0
plotheight = figwidth/Ncols
H = 1.0		#plotheight/figheight
W = 1.0 / Ncols
margin = 0.25
left = [W*margin, W*(1+margin), W*(2+margin)]
bottom = H*margin
width = W*(1-2*margin)
height = H*(1-2*margin)
piePos = 0

#########################################################
# Get gene regions (concatenate and remove duplicates)	#
#########################################################
def getGeneRegions(bedFile):
	bedFile = open(bedFile, 'r')
	regions = {}
	start = -1
	end = -1
	chrom = -1
	#if sys.argv[2] == 'contis':
	#bedFile.next()
	for line in bedFile:
		line = line.rstrip().split('\t')
		if len(line[0]) > 2:
			curChrom = line[0][3:]
		else:
			curChrom = line[0]
		curStart = int(line[1])
		curStop = int(line[2])
		associatedGeneID = line[3]
		if len(line) > 4:
			ensemblExonID = line[4]
			exon_rank = line[5]
		if associatedGeneID not in regions:
			regions[associatedGeneID] = []
		if len(line) > 4:
			regions[associatedGeneID].append([curChrom, curStart, curStop, ensemblExonID, exon_rank]) 
		else:
			regions[associatedGeneID].append([curChrom, curStart, curStop])
	sorted(regions.items(), key=lambda e: e[1])
	'''	
	else:		# No kit specified; use Ensemble coding regions database
		for line in bedFile:
			line = line.rsplit('\t')
			if line[0][2:5] == 'CHR':
				curChrom = line[0][5:]
				curChrom = curChrom.split('_')[0]
				while not curChrom.isdigit():
					curChrom = curChrom[:-1]		# Only numbered chromosomes, this excludes M, X and Y chromosomes (only M is actually missed, X and Y are not in the ensemble bed file)
				curStart = line[3]
				curStop = line[4]
				print "boe"
	'''
	return regions	

#################################################################
# Make a table of the genes that are not sufficiently covered	#
#################################################################
def LowCoveredExonsTable(depthIn, bedIn, generalBed, depth_threshold, panel):
	generalBed = open(generalBed, 'r')
	general = {}	
	if panel:
		listCoreB = []
		dataCoreB = []
		coreBFile = open(panel + "_coreB.txt", 'r')
		for line in coreBFile:
			listCoreB.append(line.rstrip())
	for line in generalBed:
		line = line.rstrip().split('\t')
		general[line[3]] = [line[0],int(line[1]),int(line[2]),line[4], line[5]]	
	lowCoveredExons = {}
	nonCoveredExons = {}
	noOfGenes = 0
	figList = []
	for key in sorted(bedIn):		#Voor alle genen
		curGeneCov = []
		for value in bedIn[key]:	#Voor alle exonen in dit gen
			#print value
			curExonCov = []
			chrom = value[0]
			start_exon = int(value[1])
			stop_exon = int(value[2])
			rank = value[4]
			if start_exon > stop_exon:		# For one exon the start and stop seems to be inverted, switch<: LELIJK
				temp = stop_exon
				stop_exon = start_exon
				start_exon = temp
			for entry in depthIn[chrom]:
				if entry[0] >= start_exon and entry[0] <= stop_exon:
					curExonCov.append([entry[0],entry[1],rank])
				elif entry[0] > stop_exon:
					break

			if curExonCov:
				for line in curExonCov:
					if line[1] < depth_threshold:
						if key not in lowCoveredExons:
							lowCoveredExons[key] = []
						lowCoveredExons[key].append(rank)	
						print "Coverage below threshold for gene " + str(key) + " exon " + str(rank)
						break	
			else:
				print "No data found for gene " + str(key) + " exon " + str(rank)
				if key not in nonCoveredExons:
					nonCoveredExons[key] = []
				nonCoveredExons[key].append(rank)
				i = start_exon
				while i <= stop_exon:
					curExonCov.append([i, 0, rank])
					i = i + 1
			curGeneCov.append(curExonCov)
	
		if key not in general:
			print 'Skipping: ' + str(key)
			continue
			
		if key in listCoreB:	#Calculate core B genelist horizontal coverage
			totalCount = 0
			countOk30 = 0
			countOk20 = 0
			countOk10 = 0
			thresholdHorizontalCov = 30
			for exone in curGeneCov:
				totalCount = totalCount + len(exone)
				for nucleotides in exone:
					if int(nucleotides[1]) >= thresholdHorizontalCov:
						countOk30 += 1
					if int(nucleotides[1]) > (thresholdHorizontalCov - 10):
						countOk20 += 1
					if int(nucleotides[1]) > (thresholdHorizontalCov - 20):
						countOk10 += 1
			dataCoreB.append([str(key), str(round(float(countOk30)/totalCount*100,2)) + " %", str(round(float(countOk20)/totalCount*100,2)) + " %", str(round(float(countOk10)/totalCount*100,2)) + " %"])

		## Plot the depth per exon
		noOfGenes = noOfGenes + 1
		xmin = general[key][1]
		xmax = general[key][2]
		plotting = []
		for depths in depthIn[chrom]:
			if depths[0] < xmin:
				continue
			elif depths[0] >= xmin and depths[0] <= xmax:
				plotting.append(depths)
			else:
				break
		fig = plt.figure(figsize=(7.5, 8.8))
		ax = fig.add_subplot(211)
		axis = []
		# Gene lies on other strand, plot in other direction...
		if general[key][3] == "-1":	
			ax.plot([-x[0] for x in plotting], [y[1] for y in plotting], '.', markersize = 2.5, color='black')
			for exon in curGeneCov:
				ax.plot([-x[0] for x in exon], [y[1] for y in exon], '.', markersize = 2.5, color = 'red')
			for value in bedIn[key]:
				plt.axvspan(int(-value[1]), int(-value[2]), facecolor = 'lightgrey', edgecolor = 'none')
				axis.append([np.mean([-value[1], -value[2]]), value[4]])	
		# Gene lies on correct strand, plot...
		else:
			ex = ax.plot([x[0] for x in plotting], [y[1] for y in plotting], '.', markersize = 2.5, color='black')
			for exon in curGeneCov:
				ax.plot([x[0] for x in exon], [y[1] for y in exon], '.', markersize = 2.5, color = 'red')
			for value in bedIn[key]:
				plt.axvspan(int(value[1]), int(value[2]), facecolor = 'lightgrey', edgecolor = 'none')
				axis.append([np.mean([value[1], value[2]]), value[4]])
		ax.axhline(depth_threshold, color = 'b', linestyle='dashed')
		plt.xticks([x[0] for x in axis], map(str,[y[1] for y in axis]), fontsize=4)
		plt.xlabel("Exon number")
		plt.yticks(fontsize=7)
		plt.ylabel("Depth")
		#plt.ylim(0, 200)
		plt.title(str(key) + " - " + str(general[key][4]), fontsize=11)
		
		#Plot information for each seperate exon
		ax = fig.add_subplot(212)
		plotAllExons = []
		for exon in curGeneCov:
			plotAllExons.append([exon[0][2], np.median([x[1] for x in exon]), np.std([x[1] for x in exon]), min([x[1] for x in exon])])
		arr = np.arange(0, len(plotAllExons))
		if plotAllExons[0][0] > plotAllExons[-1][0]:
			arr = arr[::-1]
		ax.bar(arr, [y[1] for y in plotAllExons], yerr=[ye[2] for ye in plotAllExons], capsize=2, color='lightgrey', ecolor='black', align = 'center')
		i = 0
		for minima in [x[3] for x in plotAllExons]:
			line_bar = [[arr[i]-0.4,minima],[arr[i]+0.4,minima]]
			ax.plot([x[0] for x in line_bar],[y[1] for y in line_bar], color = 'r')
			i = i+1
		ax.axhline(depth_threshold, color = 'r', linestyle='dashed')
		plt.xticks(arr, [x[0] for x in plotAllExons], fontsize = 7, rotation=90)
		plt.xlim(min(arr)-1, max(arr)+1)
		plt.ylim(ymin=0)
		plt.xlabel("Exon number")
		plt.yticks(fontsize=7)
		plt.ylabel("Median Depth")
		
		savefig(f.pathSample + '_' + str(noOfGenes) + '.png', bbox_inches='tight', dpi=300)
		figList.append(f.pathSample + '_' + str(noOfGenes) + '.png')
	
	# Here the data for the low-coverage table is filled in and printed. 
	styleSheet = getSampleStyleSheet()
	styleB = styleSheet["BodyText"]
	length = len(lowCoveredExons) + 1		
	data = []
	dataLine = None
	height_rowsl = 0.2
	for key in sorted(lowCoveredExons):
		printExonString = ""
		if len(lowCoveredExons[key]) > 13 and height_rowsl == 0.2:
			height_rowsl = 0.35
		if len(lowCoveredExons[key]) > 30:
			height_rowsl = 0.5
		for l in sorted(lowCoveredExons[key]):
			printExonString = printExonString + "  " + str(l)
		dataLine = [str(key), Paragraph(printExonString, styleB)]
		data.append(dataLine)
	if data:
		data.insert(0, [Paragraph('''<b>''' + "List of exons covered less than " + str(depth_threshold) + " x" + '''</b>''', styleB), Paragraph('''<b>''' + "All genes" + '''</b>''', styleB)])
		t = Table(data, 2 * [3 * inch], length * [height_rowsl * inch])
		t.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
						('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
						('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))	
							
	length = len(nonCoveredExons) + 1		
	data = []
	dataLine = None
	height_rowsn = 0.2
	for key in sorted(nonCoveredExons):
		printExonString = ""
		if len(nonCoveredExons[key]) > 13 and height_rowsn == 0.2:
			height_rowsn = 0.35
		if len(nonCoveredExons[key]) > 30:
			heigth_rowsn = 0.5
		for l in sorted(nonCoveredExons[key]):
			printExonString = printExonString + "  " + str(l)
		dataLine = [str(key), Paragraph(printExonString, styleB)]
		data.append(dataLine)
	if data:
		#print height_rowsn
		data.insert(0, [Paragraph('''<b>''' + "List of exons with no data" + '''</b>''', styleB), Paragraph('''<b>''' + "All genes" + '''</b>''', styleB)])
		t2 = Table(data, 2 * [3 * inch], length * [height_rowsn * inch])
		t2.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
						('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
						('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))		
		
	tpanel = []
	if panel:
		panelFile = open(panel + ".txt", 'r')
		sangerFile = open(panel + "Sanger.txt", 'r')
		listPanel = {}
		listSanger = {}	
		curPanel = ""
		for line in panelFile:
			if line[0] == "#":
				line = line.rstrip()[1:]
				listPanel[line] = []
				curPanel = line
			else:
				listPanel[curPanel].append(line.rstrip())
		#print listPanel
		#print "\n\n"
		curPanel = ""
		for line in sangerFile:
			if line[0] == "#":
				line = line.rstrip()[1:]
				listSanger[line] = []
				curPanel = line
			else:
				listSanger[curPanel].append(line.rstrip().split('\t'))
		#print listSanger
		#print "\n\n"

		data = []
		for key in listPanel:
			smallData = []
			noCov = []
			smallData.append([Paragraph('''<b>''' + str(key) + '''</b>''', styleB), Paragraph('''<b>List of exons covered less than ''' + str(depth_threshold) + ''' x</b>''', styleB)])
			for gene in listPanel[key]:
				if gene in lowCoveredExons:
					printString = ""
					for exon in lowCoveredExons[gene]:
						if gene in [x[0] for x in listSanger[key]]:
							for exones in listSanger[key]:
								if exones[0] == gene and exon in exones:
									printString = printString + '''<font color=black>''' + str(exon) + '''  </font>'''
								elif exones[0] == gene:
									printString = printString + '''<b><font color=red>''' + str(exon) + '''  </font></b>'''			
						else:
							printString = printString + '''<b><font color=red>''' + str(exon) + '''  </font></b>'''	
					printString = Paragraph(printString, styleB)
					dataLine = [str(gene), printString]
					smallData.append(dataLine)
				if gene in nonCoveredExons:
					printString = ""
					for exon in nonCoveredExons[gene]:
						if gene in [x[0] for x in listSanger[key]]:
							for exones in listSanger[key]:
								if exones[0] == gene and exon in exones:
									printString = printString + '''<font color=black>''' + str(exon) + '''  </font>'''
								elif exones[0] == gene:
									printString = printString + '''<b><font color=red>''' + str(exon) + '''  </font></b>'''
									
						else:
							printString = printString + '''<b><font color=red>''' + str(exon) + '''  </font></b>'''	
					printString = Paragraph(printString, styleB)	
					dataLine = [str(gene), printString]
					noCov.append(dataLine)
			if noCov:
				smallData.append([Paragraph('''<b>''' + str(key) + '''</b>''', styleB), Paragraph('''<b>List of exons without data</b>''', styleB)])
				smallData.extend(noCov)
			#smallData.append(["",""])
			data.append(smallData)
		if data:
			for datasets in data:
				length = len(datasets)
				t3 = Table(datasets, 2 * [3 * inch], length * [max(height_rowsn, height_rowsl) * inch])
				t3.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
						('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
						('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))	
				tpanel.append(t3)	
		if dataCoreB:
			dataCoreB.insert(0, [Paragraph('''<b>''' + "Core B Genelist" + '''</b>''', styleB), Paragraph('''<b>''' + "Horiz. cov. (>30x)" + '''</b>''', styleB), Paragraph('''<b>''' + "Horiz. cov. (>20x)" + '''</b>''', styleB), Paragraph('''<b>''' + "Horiz. cov. (>10x)" + '''</b>''', styleB)])
			length = len(dataCoreB)
			t4 = Table(dataCoreB, 4 * [1.5 * inch], length * [0.2 * inch])
			t4.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
						('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
						('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))			
			tpanel.append(t4)
					
	return t, t2, tpanel, figList

#################		
# Get depthlist #
#################
def getDepthList(depthFile):
	depthIn = open(depthFile, 'r')
	depthList = {}
	for line in depthIn:
		line = line.rstrip().rsplit('\t')
		chrom = line[0][3:]
		pos = int(line[1])
		depth = int(line[2])
		if chrom not in depthList:
			depthList[chrom] = []
		depthList[chrom].append([pos,depth])
	return depthList

	
#################################################
# Generate cumulative array from depth files	#
#################################################
def processDepth(depthList, bedRegions):
	array = []
	curDepth = {}
	'''
	# Make array not targeted
	for line in depthIn:
		line = line.rstrip().rsplit('\t')
		chrom = line[0][3:]
		pos = int(line[1])
		depth = int(line[2])
		if depth >= len(array):
			array.extend([0]*(depth - len(array) + 1))
		array[depth] += 1
	'''
	# Make array only targeted
	for gene in sorted(bedRegions):		#Voor alle genen
		for exon in bedRegions[gene]:	#Voor alle exonen in dit gen
			chrom = exon[0]
			start_exon = int(exon[1])
			stop_exon = int(exon[2])
			#rank = exon[4]
			for entry in depthList[chrom]:
				if entry[0] >= start_exon and entry[0] <= stop_exon:
					pos = entry[0]
					depth = entry[1]
					if depth >= len(array):
						array.extend([0]*(depth - len(array) + 1))
					array[depth] += 1
					if chrom not in curDepth:
						curDepth[chrom] = []
					curDepth[chrom].append([pos, depth])
					
				elif entry[0] > stop_exon:
					break
	'''
	meanMed = []
	for chrom in depthList.keys():
		for base in curDepth[chrom]:
			meanMed.append(base[-1])
	mean = np.mean(meanMed)
	median = np.median(meanMed)	
	print "old mean " + str(mean)
	print "old median " + str(median)
	
	plotDepths(array, mean, median)
	'''
	plotDepths(array)
	
#############################################################
# Generate Distribution of Depth and Cumulative Depth plots	#
#############################################################
	
def plotDepths(array):
	cumulArray = []
	maxCumul = None
	for x in range(len(array)):
		if not maxCumul:
			maxCumul = sum(array[x:])
		cumulArray.append(float(sum(array[x:])) / maxCumul * 100)
	for x in range(len(cumulArray)):
		if cumulArray[x] < 0.1:
			maxX = x - 1
			maxX = int(math.ceil(maxX / 100.0)) * 100
			break
	
	#print maxX
	#if maxX > 800:
	#	maxX = 800	
		
	sumt = 0
	amount = 0
	for i, num in enumerate(array):
		sumt += num * i
		amount += num
	
	mean = sumt / float(amount)
	if amount & 1:
		findMed = (amount + 1) / 2	#odd amount
	else:
		findMed = np.mean([(float(amount) / 2), ((float(amount) + 1) / 2)])	#even amount
	current = 0
	median = None
	for i, num in enumerate(array):
		if current < findMed:
			current += num
		else:
			median = i-1
			break
	
	array = array[1:maxX+1]

	fig = plt.figure(figsize=(8.0, 5.0))
	ax = fig.add_subplot(211)
	ax.axvline(median+0.5, color = 'r', label = 'Median')
	median = int(median)
	
	if median > maxX:
		annotate('Median: ' + str(median), xy=(median+0.5, maxX), xytext=(median+4, maxX), color = 'r')  
	else:
		annotate('Median: ' + str(median), xy=(median+0.5, array[median]), xytext=(median+4, array[median]), color = 'r')  
	ax.axvline(mean+0.5, color = '#E225D1', label = 'Mean')
	if mean > maxX:
		annotate('Mean: ' + str(mean), xy=(mean+0.5, maxX), xytext=(mean+4, maxX*0.9), color = '#E225D1') 
	else:
		annotate('Mean: ' + str(mean), xy=(mean+0.5, array[int(mean)]), xytext=(mean+4, array[int(mean)]*0.9), color = '#E225D1') 
	totalLength = len(array) 
	nrBins = 100
	if len(array) < nrBins:
		nrBins = len(array)
	itemsPerBin = (totalLength / nrBins)
	arrayBinned = []
	for i in range(nrBins):
		arrayBinned.append(np.mean(array[i * itemsPerBin : (i + 1) * itemsPerBin]))
	#ind = np.arange(0, len(arrayBinned)*itemsPerBin, itemsPerBin)
	ind = np.arange(1, len(arrayBinned)*itemsPerBin+1, itemsPerBin)
	if ind[-1] < maxX:
		width = maxX / nrBins
	else:
		width = maxX / nonzero(ind == maxX)[0][0]
	rects1 = ax.bar(ind, arrayBinned, width, color='grey')
	ax.set_ylabel('Number of bases')
	ax.set_xlabel('Depth')
	ax.set_title('Distribution of Depth (On target)')
	ax.set_xlim(0, maxX)
	ax.locator_params(nbins=6)

	ind = np.arange(len(cumulArray))  # the x locations for the groups
	ax = fig.add_subplot(212)
	rects1 = ax.plot(ind, cumulArray)
	ax.set_ylabel('% of cumulative depth')
	ax.set_xlabel('Depth')
	ax.set_title('Cumulative Depth (On target)')
	ax.set_xlim(0, maxX)
	ax.set_ylim(0, 100)
	ax.locator_params(nbins=6)
	cov = [1, minimum_depth, aimed_depth]
	for x in cov:
		ax.vlines(x, 0, cumulArray[x], linestyles='dotted')
		ax.hlines(cumulArray[x], 0, x, linestyles='dotted')
		annotate(str(round(cumulArray[x],2)) + "% at " + str(x) + " fold", xy = (x, cumulArray[x]), xytext = (x+50, cumulArray[x]-(x/2)), fontsize='x-small')
	fig.subplots_adjust(hspace=.5)

	savefig(f.pathSample + '_depth.png', dpi=300)
	

#################################################################
# Get statistics from samtools flagstat and snps.metrics file	#
#################################################################
flagstat = open(path_to_flagstat, 'r')
flTotal = int(flagstat.next().split(" ")[0])
flDups = int(flagstat.next().split(" ")[0])
flMapped = int(flagstat.next().split(" ")[0])
flPairedSeq = int(flagstat.next().split(" ")[0])
flR1 = int(flagstat.next().split(" ")[0])
flR2 = int(flagstat.next().split(" ")[0])
flPropPair = int(flagstat.next().split(" ")[0])
flBothMapped = int(flagstat.next().split(" ")[0])
flSingle = int(flagstat.next().split(" ")[0])
flDifChr = int(flagstat.next().split(" ")[0])
flDifChrQ = int(flagstat.next().split(" ")[0])
flagstat.close()

if path_to_snps_metric:
	snps_metric = open(path_to_snps_metric, 'r') 
	lList = []
	for line in snps_metric:
		lList.append(line.rstrip().rsplit()[-1])
	snps_metric.close()
		
	visitedBases = lList[0]
	callableBases = lList[1]
	confCalledBases = lList[2]
	perc1 = lList[3]
	perc2 = lList[4]
	perc3 = lList[5]
	calledBases = lList[6]

if (generalBed and bedIn) or f.hist:
	target_flagstat = open(path_to_target_flagstat, 'r')
	target_flagstat.next()
	tDup = int(target_flagstat.next().split(" ")[0])
	tMapped = int(target_flagstat.next().split(" ")[0])
	for i in xrange(4):
		target_flagstat.next()
	tPairedMap = int(target_flagstat.next().split(" ")[0])
	tSingleMap = int(target_flagstat.next().split(" ")[0])
	target_flagstat.close()

configFile = open(path_to_config[0], 'r')
rodaVersion = configFile.next()
configFile.close()
rodaVersion = rodaVersion.split()[2:]

#########################################
# Create Mapped vs. Unmapped pie chart	#
#########################################
totalMapUnmap	=	float(flTotal)
percUnmapped	=	100 * (flTotal - flMapped) / totalMapUnmap
percUnpairMap	=	100 * flSingle / totalMapUnmap
percMappedPairs	=	100 * flBothMapped / totalMapUnmap
#unknown 		=	max(0,100 - sum([percMappedPairs,percUnpairMap,percUnmapped]))
#print sum([percMappedPairs,percUnpairMap,percUnmapped])
labels=['Mapped pairs','Mapped Singles','Unmapped']
fracs = [percMappedPairs,percUnpairMap,percUnmapped]
legend1 = [str(fracs[i])[:5] + '% ' + labels[i] for i in range(len(labels))]
explode=(0.0, 0.0, 0.0)
colors=('#34C934','#E26A06','#FF0000')
axes([left[piePos], bottom, width, height])
patches = pie(fracs, labels=labels,colors=colors, autopct='%1.f%%', shadow=False) #,  explode=explode,
title('Mapped vs. Unmapped')
legend(legend1,loc=(-0.4,-0.4))
piePos += 1

#################################
# Create Duplicates pie chart	#
#################################
totalOptDups 	=	flDups
totalDupsUndups = 	float(flMapped)		#= dus 1
totalNonDups 	= 	totalDupsUndups - totalOptDups
percOptDups		=	100 * (totalOptDups / totalDupsUndups)
percNonDups		=	100 * (totalNonDups / totalDupsUndups)

labels=['Unique mapped','Duplicates']
fracs = [percNonDups,percOptDups]
legend1 = [str(fracs[i])[:5] + '% ' + labels[i] for i in range(len(labels))]
explode=(0.0, 0.0)
colors=('#34C934','#FF0000')
axes([left[piePos], bottom, width, height])
patches = pie(fracs, labels=labels,colors=colors, autopct='%1.f%%', shadow=False) #,  explode=explode,
title('Duplicates')
legend(legend1,loc=(-0.4,-0.4))
piePos += 1

#############################
# Create Target pie chart	#
#############################
if bedIn or f.hist:
	totalTarget 	=	float(tMapped + (flMapped - tMapped))
	#totalTarget		=	float(tPairedMap + tSingleMap + tDup)
	#perctDups		=	100 * (tDup / totalTarget)
	perctPaired		=	100 * (tPairedMap / totalTarget)
	perctSingle		=	100 * (tSingleMap / totalTarget)
	notOnTarget		= 	100 * ((flMapped - tMapped) / totalTarget)

	#labels=['Mapped Pairs','Mapped Singles', 'Duplicates']
	labels=['Pairs on target','Singles on target', 'Not on target']
	#fracs = [perctPaired,perctSingle,perctDups]
	fracs = [perctPaired,perctSingle,notOnTarget]
	legend1 = [str(fracs[i])[:5] + '% ' + labels[i] for i in range(len(labels))]
	explode=(0.0, 0.0, 0.0)
	colors=('#34C934','#E26A06','#FF0000')
	axes([left[piePos], bottom, width, height])
	patches = pie(fracs, labels=labels,colors=colors, autopct='%1.f%%', shadow=False) #,  explode=explode,
	title('Mapped reads')
	legend(legend1,loc=(-0.4,-0.4))
	piePos += 1

savefig(f.pathSample + '_pie.png')


#########################################
# Print metrics statistics to tables	#
#########################################
elements.append(Paragraph('RoDa Version ' + ' '.join(rodaVersion) + ' Quality report', styleHeading))
elements.append(Paragraph('Sample ' + sampleName, styleHeading))
elements.append(Spacer(inch, 0.25*inch)	)					
elements.append(Paragraph('Mapping quality report', styleNormal))

data = [['Total reads', 				str(flTotal), 				str(100) + ' %'],
		['Mapped paired reads', 		str(flBothMapped),			str(round((flBothMapped / float(flTotal) * 100),2)) + ' %'],
		['Mapped unpaired reads', 		str(flSingle),				str(round((flSingle / float(flTotal) * 100),2)) + ' %'],
		['Unmapped reads', 				str(flTotal - flMapped), 	str(round(((flTotal - flMapped) / float(flTotal) * 100),2)) + ' %'],
		['Total optical duplicates',	str(flDups),				str(round((flDups / float(flTotal) * 100),2)) + ' %']]
if bedIn or f.hist:
	data.append(['Total mapped on target', 		str(tMapped), 				str(round(((tMapped / totalTarget) * 100), 2)) + ' %'])
	data.append(['Not on target', 				str(flMapped - tMapped), 	str(round((((flMapped - tMapped)/totalTarget) * 100), 2)) + ' %'])
	#data.append(['Total on target', 			str(totalTarget), 				str(round(((totalTarget /flTotal) * 100), 2)) + ' %'])
	#data.append(['On target minus duplicates',	str(perctPaired + perctSingle),	str(round((((perctPaired + perctSingle)/flTotal) * 100), 2)) + ' %'])
	#data.append(['Not on target', 				str(flTotal - totalTarget), 	str(round((((flTotal - totalTarget)/flTotal) * 100), 2)) + ' %'])
t = Table(data, 3 * [2 * inch], len(data) * [0.2 * inch])
t.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
						('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
						('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))
						
						
elements.append(t)
elements.append(Image(f.pathSample + '_pie.png', 7.5*inch, 2.5*inch))

if generalBed and not bedIn and f.hist:
	generalRegions = getGeneRegions(generalBed)
	array = []
	hist = open(path_to_histFile, 'r')
	over30Cov = int(hist.next().rstrip())
	totalBasesExome = int(hist.next().rstrip())
	for line in hist:
		line = line.rstrip().split("\t")
		array.append(int(line[2]))
	#print array
	plotDepths(array)
	elements.append(Image(f.pathSample + '_depth.png', 8*inch, 5*inch))
	elements.append(Spacer(inch, 1*inch))
	elements.append(Paragraph('Total bases exome:\t' + str(totalBasesExome), styleNormal))
	elements.append(Paragraph('Total covered bases (>30x):\t' + str(over30Cov), styleNormal))
	elements.append(Paragraph('Horizontal coverage target (>30x):\t' + str(round(over30Cov/float(totalBasesExome)*100, 2)) + " %", styleNormal))
	elements.append(PageBreak())


if bedIn and generalBed:
	bedRegions = getGeneRegions(bedIn)
	depthProcessed = getDepthList(depthFile)
	processDepth(depthProcessed, bedRegions)
	elements.append(Image(f.pathSample + '_depth.png', 8*inch, 5*inch))
	elements.append(PageBreak())
	t, t2, tpanel, figList = LowCoveredExonsTable(depthProcessed, bedRegions, generalBed, depth_threshold, f.panel)
	elements.append(t)
	elements.append(t2)
	elements.append(Spacer(inch, 1*inch))
	if tpanel:
		for ts in tpanel:
			elements.append(ts)
			elements.append(Spacer(inch, 0.35*inch))
	for figure in figList:
		elements.append(Image(figure, 7.5*inch, 8.8*inch))

if path_to_snps_metric:
	data = [['Visited bases',								str(visitedBases)],
			['Callable bases',								str(callableBases)],
			['Confidently called bases',					str(confCalledBases)],
			['Total called variants',						str(calledBases)],
			['Callable bases of all loci',					str(perc1) + ' %'],
			['Confidently called bases of all loci',		str(perc2) + ' %'],
			['Confidently called bases of callable loci',	str(perc3) + ' %']]
	t = Table(data,2 * [3 * inch], 7 * [0.2 * inch])
	t.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
				('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
				('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))
	elements.append(Paragraph('SNP/indel caller stats', styleNormal))
	elements.append(t)

#############################
# Export the pdf document	#
#############################
outputPDF.build(elements)
if bedIn:
	for files in figList:
		os.remove(files)

