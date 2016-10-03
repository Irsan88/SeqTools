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
import datetime
	
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
	#path_to_snps_metric = f.pathSample + ".stat.snps.metrics"		#".snps.metrics"
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
styleNormal2 = styles["Normal"]
styleNormal2.alignment = 2
elements = []

Ncols = 3
Nrows = 2
edge = 2.5
figheight = edge * Nrows   # inches
figwidth = edge * Ncols    # inches
figure(1, figsize=(figwidth, figheight))
rcParams['font.size'] = 6.0
rcParams['axes.titlesize'] = 12.0
rcParams['xtick.labelsize'] = 7.0
rcParams['legend.fontsize'] = 5.0
plotheight = figwidth/Ncols
H = 1.0	/ Nrows	#plotheight/figheight
W = 1.0 / Ncols
margin = 0.25
left = [W*margin, W*(1+margin), W*(2+margin), W*(1+margin)]
bottom = [H*(margin+1), H*(margin+1), H*(margin+1), H*margin]
widthplot = W*(1-2*margin)
heightplot = H*(1-2*margin)

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
		if len(line) > 5:
			ensemblExonID = line[4]
			exon_rank = line[5]
		if associatedGeneID not in regions:
			regions[associatedGeneID] = []
		if len(line) > 5:
			regions[associatedGeneID].append([curChrom, curStart, curStop, ensemblExonID, exon_rank]) 
		else:
			regions[associatedGeneID].append([curChrom, curStart, curStop])
	sorted(regions.items(), key=lambda e: e[1])
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
						#print "Coverage below threshold for gene " + str(key) + " exon " + str(rank)
						break	
			else:
				#print "No data found for gene " + str(key) + " exon " + str(rank)
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
		ax.axhline(depth_threshold, color = 'b', linestyle='dashed')
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
	countLowCov = 0
	countNoCov = 0
	for key in sorted(lowCoveredExons):
		printExonString = ""
		if len(lowCoveredExons[key]) > 13 and height_rowsl == 0.2:
			height_rowsl = 0.35
		if len(lowCoveredExons[key]) > 30:
			height_rowsl = 0.5
		for l in sorted(lowCoveredExons[key]):
			countLowCov += 1
			printExonString = printExonString + "  " + str(l)
		dataLine = [str(key), Paragraph(printExonString, styleB)]
		data.append(dataLine)
	if data:
		data.insert(0, [Paragraph('''<b>''' + "ALL GENES" + '''</b>''', styleB), Paragraph('''<b>''' + "List of exons covered less than " + str(depth_threshold) + " x" + '''</b>''', styleB)])
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
			countNoCov += 1
			printExonString = printExonString + "  " + str(l)
		dataLine = [str(key), Paragraph(printExonString, styleB)]
		data.append(dataLine)
	if data:
		#print height_rowsn
		data.insert(0, [Paragraph('''<b>''' + "ALL GENES" + '''</b>''', styleB), Paragraph('''<b>''' + "List of exons with no data" + '''</b>''', styleB)])
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

		countAlwaysSanger = 0
		data = []
		countRed = 0
		for key in listPanel:
			smallData = []
			noCov = []
			alwaysSanger = []
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
									countRed += 1
						else:
							printString = printString + '''<b><font color=red>''' + str(exon) + '''  </font></b>'''
							countRed += 1
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
									countRed += 1
						else:
							printString = printString + '''<b><font color=red>''' + str(exon) + '''  </font></b>'''
							countRed += 1	
					printString = Paragraph(printString, styleB)	
					dataLine = [str(gene), printString]
					noCov.append(dataLine)
				if gene in [x[0] for x in listSanger[key]]:
					printString = ""
					for exones in listSanger[key]:
						if exones[0] == gene:
							for ex in exones:
								if ex[0:1] == "!":
									printString = printString + '''<font color=green>''' + str(ex[1:]) + '''  </font>'''
									countAlwaysSanger += 1
					if printString:
						printString = Paragraph(printString, styleB)
						dataLine = [str(gene), printString]
						alwaysSanger.append(dataLine)
					
				
							
			if noCov:
				smallData.append([Paragraph('''<b>''' + str(key) + '''</b>''', styleB), Paragraph('''<b>List of exons with no data</b>''', styleB)])
				smallData.extend(noCov)
			if alwaysSanger:
				smallData.append([Paragraph('''<b>''' + str(key) + '''</b>''', styleB), Paragraph('''<b>Exons that are always Sanger sequenced</b>''', styleB)])
				smallData.extend(alwaysSanger)

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
			dataCoreB.insert(0, [Paragraph('''<b>''' + "Core B Genelist" + '''</b>''', styleB), Paragraph('''<b>''' + "Horiz. cov. (>=30x)" + '''</b>''', styleB), Paragraph('''<b>''' + "Horiz. cov. (>=20x)" + '''</b>''', styleB), Paragraph('''<b>''' + "Horiz. cov. (>=10x)" + '''</b>''', styleB)])
			length = len(dataCoreB)
			t4 = Table(dataCoreB, 4 * [1.5 * inch], length * [0.2 * inch])
			t4.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
						('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
						('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))			
			tpanel.append(t4)
	print "Low covered exons\tNon covered exons\tAlways Sanger sequenced"
	print str(countLowCov) + "\t" + str(countNoCov) + "\t" + str(countAlwaysSanger) 	
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
# If file is available: plot MQ distribution					#
#################################################################
def plotMQDist():
	try:
		distFile = open(f.pathSample + ".stat.mqdist", 'r')
		countTotalReads = 0
		distX = []
		distY = []
		noUnder60 = 0
		total = 0
		for line in distFile:
			line = line.rstrip().split("\t")
			distX.append(int(line[0]))
			distY.append(int(line[1]))
			countTotalReads += int(line[1])
			if int(line[0]) < 60:
				noUnder60 += int(line[1])
			total += int(line[1])
		fig = plt.figure(figsize=(7.5, 3))
		ax = fig.add_subplot(111)
		ax.bar(distX, distY)
		ax.set_ylabel('Count')
		ax.set_xlabel('Mapping Quality')
		ax.set_title('Distribution of Mapping Quality')
		ax.set_xlim(0,max(distX))
		percentage = str(round(noUnder60 / float(total) * 100, 2))
		ax.hlines(noUnder60, 0, max(distX), linestyles='dotted', label=percentage)
		ax.text(20, noUnder60 *1.1, r'% MQ < 60: ' + percentage + '%', fontsize=10)
		savefig(f.pathSample + '_mqdist.png', dpi=300)
		#print countTotalReads
		return True
	except IOError:
		print "No MQ dist file found."
		return False

#################################################################
# For exomes: plot distribution hom/het calls					#
#################################################################
def plotHomHet():
	try:
		hethomFile = open(f.pathSample + ".snps.filt.RD.vcf", 'r')
		positions = {}
		allChroms = []
		lengthContigs = {}
		longHomRegs = {}

		for line in hethomFile:
			if line[0] == '#':
				if line[0:8] == '##contig':
					line = line.rstrip().split(',')
					chrom = line[0].split('=')[-1]
					length = line[1].split('=')[-1]
					lengthContigs[chrom] = length
				continue
			line = line.rstrip().split('\t')
			if line[6] != 'PASS':
				continue
			if line[0] not in positions:
				positions[line[0]] = []
				allChroms.append(line[0])
			alleleFreq = line[7].split(';')[1]
			alleleFreq = alleleFreq.split('=')
			if alleleFreq[0] != 'AF':
				print "Warning: AF not at expected position: \t" + str(alleleFreq)
			positions[line[0]].append([line[1], alleleFreq[1]])
			
		subplots = [0, 411, 412, 413, 414]
		loop = 0
		nrFig = 0
		allChroms.sort(key=lambda x: x[0])
		fig = plt.figure(figsize=(20, 10), dpi=300)
		
		for chrom in allChroms:
			curROH = None
			loop = loop + 1
			if loop == 5:
				nrFig = nrFig + 1
				plt.savefig(f.pathSample + "_" + str(nrFig) + '.jpg')
				fig = plt.figure(figsize=(20, 10), dpi=300)
				loop = 1
			posx = []
			posy = []
			for pos in positions[chrom]:
				freqs = pos[1].split(',')
				lengthfreqs = len(freqs)
				if float(freqs[0]) > 0.75:
					#print "HOM!"
					if not curROH:
						#print 'making curROH'
						curROH = [int(pos[0]),None,0]
						#print 'filling in second pos ROH'
					curROH[1] = int(pos[0])
					curROH[2] += 1 
				else:
					if curROH:
						#print 'Found Het, break off curROH'
						#if curROH[2] > 10 or ((curROH[1] - curROH[0]) > 1000000 and curROH[2] > 25):
						if (curROH[1] - curROH[0]) > 500000 and curROH[2] > 20:
							#print 'curROH is big!'
							if chrom not in longHomRegs:
								longHomRegs[chrom] = []
							longHomRegs[chrom].append(curROH)
						curROH = None
				
				for i in range(0,lengthfreqs):
					posx.append(int(pos[0]))
					posy.append(float(freqs[i]))	
			ax = fig.add_subplot(subplots[loop])
			if chrom in longHomRegs:
				print chrom
				print longHomRegs[chrom]
				for tups in longHomRegs[chrom]:
					if int(tups[2]) > 25:
						plt.axvspan(int(tups[0]), int(tups[1]), facecolor = 'red', alpha=0.2, edgecolor = 'none')
					else:
						plt.axvspan(int(tups[0]), int(tups[1]), facecolor = 'blue', alpha=0.2, edgecolor = 'none')
			ax.scatter(posx, posy, s=2, facecolor='0.5', lw = 0)
			ax.set_ylabel('AF ' + str(chrom))
			ax.set_ylim(-0.1,1.1)
			ax.set_xlim(0, int(lengthContigs[chrom]))
			locs,labels = xticks()
			xticks(locs, map(lambda posx: "%0.0f" % posx, locs))
			for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
				item.set_fontsize(10)
		nrFig = nrFig + 1
		plt.savefig(f.pathSample + "_" + str(nrFig) + '.jpg')
		return True
	
	except IOError:
		print "No MQ dist file found."
		return False	

#################################################################
# Get statistics from samtools flagstat and snps.metrics file	#
#################################################################
flagstat = open(path_to_flagstat, 'r')
flTotal = float(flagstat.next().split(" ")[0])
flDups = int(flagstat.next().split(" ")[0])
flMapped = float(flagstat.next().split(" ")[0])
flPairedSeq = int(flagstat.next().split(" ")[0])
flR1 = int(flagstat.next().split(" ")[0])
flR2 = int(flagstat.next().split(" ")[0])
flPropPair = int(flagstat.next().split(" ")[0])
flBothMapped = int(flagstat.next().split(" ")[0])
flSingle = int(flagstat.next().split(" ")[0])
flDifChr = int(flagstat.next().split(" ")[0])
flDifChrQ = int(flagstat.next().split(" ")[0])
flagstat.next()
flMQ0 = int(flagstat.next().split(" ")[-1])
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

target_flagstat = open(path_to_target_flagstat, 'r')
tTotal = int(target_flagstat.next().split(" ")[0])
tDups = int(target_flagstat.next().split(" ")[0])
tMapped = int(target_flagstat.next().split(" ")[0])
tPairedSeq = int(target_flagstat.next().split(" ")[0])
tR1 = int(target_flagstat.next().split(" ")[0])
tR2 = int(target_flagstat.next().split(" ")[0])
tPropPair = int(target_flagstat.next().split(" ")[0])
tBothMapped = int(target_flagstat.next().split(" ")[0])
tSingle = int(target_flagstat.next().split(" ")[0])
tDifChr = int(target_flagstat.next().split(" ")[0])
tDifChrQ = int(target_flagstat.next().split(" ")[0])
tMappedOnRS = int(target_flagstat.next().split(" ")[-1])
tMappedOnFS = int(target_flagstat.next().split(" ")[-1])
target_flagstat.close()

configFile = open(path_to_config[0], 'r')
rodaVersion = configFile.next().split()[2:]
dateRun = configFile.next().split()[-1]
configFile.close()

#########################################
# Create Mapped vs. Unmapped pie chart	#
#########################################
flUnmapped		= 	flTotal - flMapped
percUnmapped	=	100 * flUnmapped / flTotal
percSingleton	=	100 * flSingle / flTotal
percMappedPairs	=	100 * flBothMapped / flTotal

labels=['Mapped pairs','Mapped Singles','Unmapped']
fracs = [percMappedPairs,percSingleton,percUnmapped]
legend1 = [str(fracs[i])[:5] + '% ' + labels[i] for i in range(len(labels))]
explode=(0.0, 0.0, 0.0)
colors=('#34C934','#0000FF','#FF0000')
axes([left[piePos], bottom[piePos], widthplot, heightplot])
patches = pie(fracs, labels=labels,colors=colors, autopct='%1.f%%', shadow=False) #,  explode=explode,
title('Mapped vs. Unmapped')
legend(legend1,loc=(-0.4,-0.4))
piePos += 1

#################################
# Create Duplicates pie chart	#
#################################
totalNonDups 	= 	flMapped - flDups
percOptDups		=	100 * (flDups / flMapped)
percNonDups		=	100 * (totalNonDups / flMapped)

labels=['Unique mapped','Duplicates']
fracs = [percNonDups,percOptDups]
legend1 = [str(fracs[i])[:5] + '% ' + labels[i] for i in range(len(labels))]
explode=(0.0, 0.0)
colors=('#34C934','#FF0000')
axes([left[piePos], bottom[piePos], widthplot, heightplot])
patches = pie(fracs, labels=labels,colors=colors, autopct='%1.f%%', shadow=False) #,  explode=explode,
title('Duplicates')
legend(legend1,loc=(-0.4,-0.4))
piePos += 1

#############################
# Create Target pie chart	#
#############################
#print flMapped
#print (tSingle + (flMapped - tMapped) + tMappedOnFS + tMappedOnRS)

perctSingle		=	100 * (tSingle / flMapped)
notOnTarget		= 	100 * ((flMapped - tMapped) / flMapped)
perctMapped		=	100 * (tMapped / flMapped)

perctPairedFS	=	100 * (tMappedOnFS / flMapped)
perctPairedRS	=	100 * (tMappedOnRS / flMapped)
labels=['Pairs on target FS','Pairs on target RS','Singles on target','Not on target']
fracs = [perctPairedFS,perctPairedRS,perctSingle,notOnTarget]
legend1 = [str(fracs[i])[:5] + '% ' + labels[i] for i in range(len(labels))]
explode=(0.0, 0.0, 0.0, 0.0)
colors=('#34C934','#0000FF','#FF8000','#FF0000')
axes([left[piePos], bottom[piePos], widthplot, heightplot])
patches = pie(fracs, labels=labels,colors=colors, autopct='%1.f%%', shadow=False) #,  explode=explode,
title('Mapped reads')
legend(legend1,loc=(-0.4,-0.4))
piePos += 1

	
#############################
# Create Summary pie chart	#
#############################
sumMQOffT	 	=	100 * flMQ0 / flTotal
totalMap		= 	(flBothMapped + flSingle - flDups - flMQ0)
unmapped		=	(flTotal - (totalMap + flDups + flMQ0))
	
onTargetMap 	= 	(tBothMapped + tSingle - tDups)
mappedOffT		=	(totalMap - onTargetMap)
offTDups		= 	(flDups - tDups)
sumMappedTarget =	100 * onTargetMap / flTotal
sumMappedOffT	=	100 * mappedOffT / flTotal
sumUnmapped		=	100 * unmapped / flTotal
sumDupsOffT		=	100 * (flDups - tDups) / flTotal
sumDupsTarget	=	100 * tDups / flTotal

labels=['Unmapped','Duplicates on target','Duplicates off target','Mapped on target','Mapped off target','MQ0']
fracs = [sumUnmapped,sumDupsTarget,sumDupsOffT,sumMappedTarget,sumMappedOffT,sumMQOffT]
#fracs = [percUnmapped,percOptDups,percOnTargetMap,percMQ0,percOffTarget]
legend1 = [str(fracs[i])[:5] + '% ' + labels[i] for i in range(len(labels))]
explode=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
colors=('#FF0000','#0B610B','#FFFF00','#34C934','#FFA500','#FACC2E')
axes([left[piePos], bottom[piePos], widthplot*1.2, heightplot*1.2])
patches = pie(fracs, labels=labels,colors=colors, autopct='%1.f%%', shadow=False) #,  explode=explode,
title('Summary')
legend(legend1,loc=(-0.4,-0.4))
piePos += 1

savefig(f.pathSample + '_pie.png', dpi=300)

print "Total Reads\tUnmapped\tDuplicates on target\tDuplicates off target\tMapped on target\tMapped off target\tMQ0"
print str(int(flTotal)) + '\t' + str(int(unmapped)) + '\t' + str(tDups) + '\t' + str(flDups - tDups) + '\t' + str(onTargetMap) + '\t' + str(mappedOffT) + '\t' + str(flMQ0)

#########################################
# Print metrics statistics to tables	#
#########################################
elements.append(Paragraph('RoDa Version ' + ' '.join(rodaVersion) + ' Quality report', styleHeading))
elements.append(Paragraph('Sample ' + sampleName, styleHeading))
elements.append(Spacer(inch, 0.1*inch)	)
now = datetime.datetime.now()
elements.append(Paragraph('Date analysis performed:\t' + str(dateRun), styleNormal2))
elements.append(Paragraph('Date printing report:\t' + now.strftime("%d/%m/%y"), styleNormal2))
elements.append(Spacer(inch, 0.25*inch)	)	
elements.append(Paragraph('Mapping quality report', styleNormal))				

data = [['Total reads', 				str(int(flTotal)), 				str(100) + ' %'],
		['Mapped paired reads', 		str(flBothMapped),			str(round((flBothMapped / float(flTotal) * 100),2)) + ' %'],
		['Mapped unpaired reads', 		str(flSingle),				str(round((flSingle / float(flTotal) * 100),2)) + ' %'],
		['Unmapped reads', 				str(int(unmapped)), 				str(round((unmapped / float(flTotal) * 100),2)) + ' %'],
		['Total optical duplicates',	str(flDups),				str(round((flDups / float(flTotal) * 100),2)) + ' %'],
		['Total mapped on target', 		str(onTargetMap), 			str(round(((onTargetMap / flTotal) * 100), 2)) + ' %'],
		['Not on target', 				str(mappedOffT), 			str(round(((mappedOffT/flTotal) * 100), 2)) + ' %']]
	
t = Table(data, 3 * [2 * inch], len(data) * [0.2 * inch])
t.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
						('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
						('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))
						
						
elements.append(t)
elements.append(Image(f.pathSample + '_pie.png', figwidth*inch, figheight*inch))

rcParams['font.size'] = 9.0
rcParams['xtick.labelsize'] = 9.0

if f.hist:
	#generalRegions = getGeneRegions(generalBed)
	OLD = False
	generalRegions = getGeneRegions(bedIn)
	array = []
	hist = open(path_to_histFile, 'r')
	over1Cov = int(hist.next().rstrip())
	over10Cov = int(hist.next().rstrip())
	line = hist.next().rstrip()
	if line[0:3] == "all":
		#old version of hist
		print "Old version of hist detected!"
		OLD = True
		line = line.split("\t")
		array.append(int(line[2]))
		over30Cov = over1Cov
		totalBasesExome = over10Cov
	else:
		#new version of hist
		over30Cov = int(line)
		over50Cov = int(hist.next().rstrip())
		totalBasesExome = int(hist.next().rstrip())
	for line in hist:
		line = line.rstrip().split("\t")
		array.append(int(line[2]))
	plotDepths(array)
	elements.append(Image(f.pathSample + '_depth.png', 8*inch, 5*inch))
	elements.append(Spacer(inch, 1*inch))
	if not OLD:
		data = [[Paragraph('<b>Description</b>', styleNormal),	Paragraph('<b>Total number of bases</b>', styleNormal),			Paragraph('<b>Horizontal coverage</b>', styleNormal)],
			['Whole exome', 	str(totalBasesExome),	'-'],
			['>=1x, MQ>=10', 		str(over1Cov),			str(round(over1Cov/float(totalBasesExome)*100, 2)) + " %"],	
			['>=10x, MQ>=10', 	str(over10Cov),			str(round(over10Cov/float(totalBasesExome)*100, 2)) + " %"],
			['>=30x, MQ>=10', 	str(over30Cov),			str(round(over30Cov/float(totalBasesExome)*100, 2)) + " %"],
			['>=50x, MQ>=10', 	str(over50Cov),			str(round(over50Cov/float(totalBasesExome)*100, 2)) + " %"]]
		t = Table(data,3 * [2 * inch], 6 * [0.2 * inch])
		t.setStyle(TableStyle([	('INNERGRID', (0,0), (-1,-1), 0.25, ncolors.black),
				('BOX', (0,0), (-1,-1), 0.25, ncolors.black),
				('VALIGN', (0,0), (-1,-1), 'MIDDLE')]))
		elements.append(Paragraph('Coverage statistics', styleNormal))
		elements.append(t)
	else:
		elements.append(Spacer(inch, 1*inch))
		elements.append(Paragraph('Total bases exome:\t' + str(totalBasesExome), styleNormal))
		elements.append(Paragraph('Total covered bases (>30x, MQ>10):\t' + str(over30Cov), styleNormal))
		elements.append(Paragraph('Horizontal coverage target (>30x, MQ>10):\t' + str(round(over30Cov/float(totalBasesExome)*100, 2)) + " %", styleNormal))

	
	elements.append(PageBreak())
	plotHomHet()


if not f.hist:
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
	
if plotMQDist() is True:
	elements.append(Spacer(inch, 0.35*inch))
	elements.append(Image(f.pathSample + '_mqdist.png', 7.5*inch, 3*inch))

#############################
# Export the pdf document	#
#############################
outputPDF.build(elements)
if not f.hist:
	for files in figList:
		os.remove(files)

