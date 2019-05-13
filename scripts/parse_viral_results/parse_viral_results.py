#!/usr/bin/python
import sys
import os
import glob
import re
import string
#import common.json

def replaceOldDB(filename, container): ###Function for Methaphlan2 container
	skip = 0 #Skipping header
	d = {} #Dictonary container
	with open(filename) as f:
		for line in f:
			(key, val) = line.split()
			d[str(key)] = str(val)
	f.close()
	for c in container:
		if skip != 0 and c[1] in d:
			c[1] = d[str(c[1])]
		skip = 1 # value after skipping header
	return container

def taxonLevelString(arr, select):
	strn = ""
	tmp = []
	#Select 1 for kraken and select 2 for methaphlan
	if select == 1:
		for a in arr:
			if strn == "":
				strn = strn + a
			else:
				strn = strn + "->" + a
	if select == 2:
		for a in arr:
			if "f__" in a:
				tmp.append(a.replace("f__", ""))
			elif "g__" in a:
				tmp.append(a.replace("g__", ""))
			elif "s__" in a:
				tmp.append(a.replace("s__", ""))
			elif "t__" in a:
				tmp.append(a.replace("t__", ""))
		for a in tmp:
			if strn == "":
				strn = strn + a
			else:
				strn = strn + "->" + a
	return strn
		
def CreateFile(fileName):
	name = "Combined_Results_"+fileName+".txt"
	f= open(name,"w+")
	f.close() 
	return None

def WriteToFile(resultArray, fileName):
	name = "Combined_Results_"+fileName+".txt"
	tmp = ""
	f=open(name, "a+")
	for r in resultArray:
		tmp = str(r[0])+"\t"+str(r[1])+"\t"+str(r[2])+"\t"+str(r[3])+"\t"+str(r[4])+"\t"+str(r[5])+"\t"+str(r[6])+"\n"
		f.write(tmp)
	f.close() 
	return None

def ireplace(old, new, text):
	index_l = text.lower().index(old.lower())
	return text[:index_l] + new + text[index_l + len(old):]

def parseMetaphlan2(fileInput):
	#Classification info
	"""
	Kingdom[0]: k__, Phylum[1]: p__, Class[2]: c__, Order[3]: o__, Family[4]: f__, Genus[5]: g__, Species[6]: s__, Sub_species/Strain[7]: t__
	
	Since sequence-based profiling is relative and does not provide absolute cellular abundance measures, 
	clades are hierarchically summed. Each level will sum to 100%; that is, 
	the sum of all kindom-level clades is 100%, the sum of all genus-level clades (including unclassified) is also 100%, and so forth. 
	OTU equivalents can be extracted by using only the species-level s__ clades from this file (again, making sure to include clades unclassified at this level).
	
	Source: https://bitbucket.org/biobakery/biobakery/wiki/metaphlan2
	
	Output from function:
	['Classification', 'Name', 'Taxon Order', 'Abundance']
	['Family', 'Flaviviridae', 'k__Viruses|p__Viruses_noname|c__Viruses_noname|o__Viruses_noname|f__Flaviviridae', '100.0']  
	"""
	#Header from table
	f = open(fileInput)
	line = f.readline().split()
	tmp =  []
	prevTmp = []
	prevTmpContent = []
	container = []
	classType = []
	
	#sample ID from header
	SampleID = line[1]

	#Reading file
	addClass = ['f__','g__','s__','t__']
	ClassNames = ['Family','Genus','Species','Sub_Species']
	line = f.readline().split()
	tmpline = line[0].split('|')
	prevTmp = tmpline
	prevTmpContent = line
	prevTmpContent.insert(0, prevTmp[len(prevTmp)-1])
	run = 1
	control = 1
	while run != 0:
		line = f.readline().split()
		if len(line) == 0:
			for i in range(len(addClass)):
				if addClass[i] in prevTmpContent[0]:
					prevTmpContent[0] = prevTmpContent[0].replace(addClass[i],'')
					prevTmpContent.insert(0, ClassNames[i])
			if "virus" in str(prevTmpContent[2]) or "Virus" in str(prevTmpContent[2]):
				container.append(prevTmpContent)
			break
		tmpline = line[0].split('|')
		
		for n in prevTmp:
			if any(str(it) in str(n) for it in addClass):
				control = 0
				break
		if control != 1: 
			for i in range(len(addClass)):
				if addClass[i] in prevTmpContent[0]:
					prevTmpContent[0] = prevTmpContent[0].replace(addClass[i],'')
					prevTmpContent.insert(0, ClassNames[i])
			if "virus" in str(prevTmpContent[2]) or "Virus" in str(prevTmpContent[2]):
				container.append(prevTmpContent)
		control = 1
		prevTmp = tmpline
		prevTmpContent = line
		prevTmpContent.insert(0, prevTmp[len(prevTmp)-1])
	container.insert(0, [SampleID])
	"""
	for n in container:
		print(n)
		print('\n')
	"""
	container = replaceOldDB("DictonaryDB.txt", container)
	return container

def fixFVElist(inputList):
	#print(inputList[0])
	newlist = [inputList[0]]
	inputList = inputList[1:len(inputList)]
	check = 1
	virusName = ''
	className = ''
	while check != 0:
		if ";" not in inputList[0]:
			virusName = virusName + inputList[0]
			inputList = inputList[1:len(inputList)]
			if ";" not in inputList[0]:
				virusName = virusName +"_"
		else:
			newlist.append(virusName)
			break
	while check != 0:
		if len(inputList) != 1:
			className = className + inputList[0]
			inputList = inputList[1:len(inputList)]
			if len(inputList) != 1:
				className = className +"_"
		else:
			newlist.append(className)
			newlist.append(inputList[0])
			break
	classStr = newlist[2].split(';')
	index = int(len(classStr));
	label = " "
	if index == 8:
		classStr = classStr[4:8]
		label = "Sub_Species"
	elif index == 7:
		classStr = classStr[4:7]
		label = "Species"
	elif index == 6:
		classStr = classStr[4:6]
		label = "Genus"
	elif index == 5:
		classStr = classStr[4:5]
		label = "Family"
	elif index == 4:
		classStr = classStr[4:4]
		label = "Family"
	else:
		classStr = ['delete']
		
	#Printing for debuggin index
	#print(index)
	#print(len(classStr))
	#print(classStr)
	
	#2 new colomns added index 0 gives taxon level and index 1 gives name of taxon level same as in methaphlan
	newlist.insert(0, classStr[len(classStr)-1])
	newlist.insert(0, label)
	return newlist

def parseFVE(fileInput):
	#Parsing FastViromeExplorer output data
	#Classification info
	"""
	kingdom;	phylum;	class;	order;	family;	genus;	species
	
	Header:
	['#VirusIdentifier', 'VirusName', 'kingdom;phylum;class;order;family;genus;species', 'EstimatedAbundance']
	
	unfixed
	['NC_001710.1', 'GB', 'virus', 'C/Hepatitis', 'G', 'virus,', 'complete', 'genome', 'Unclassified;Unclassified;Unclassified;Unclassified;Flaviviridae;Pegivirus;Pegivirus', 'C', '66.0']
	
	fixed output from fixFVElist
	['NC_001710.1', 'GB_virus_C/Hepatitis_G_virus,_complete_genome', 'Unclassified;Unclassified;Unclassified;Unclassified;Flaviviridae;Pegivirus;Pegivirus_C', '66.0']
	
	new fixed index 0 = taxon level and index 1 = name of taxon level
	['Species', 'Human_mastadenovirus_C', 'AC_000008.1', 'Human_adenovirus_5,_complete_genome', 'Unclassified;Unclassified;Unclassified;Unclassified;Adenoviridae;Mastadenovirus;Human_mastadenovirus_C', '1248.75']  
	"""
	#Skipping Header
	f = open(fileInput)
	line = f.readline().split()
	#Reading file FIXA loop för att läsa filen från FVE
	container = []
	run = 1
	while run != 0:
		line = f.readline().split()
		#print(line)
		if len(line) != 0:
			line = fixFVElist(line)
			container.append(line)
		else:
			break
	return container

def FVEaddSubSpecies(prevCont):
	container = []
	tmp = []
	tmpRow = []
	tmpRow2 = []
	subSpecies = []
	check = 0
	sumValues = 0.0
	for n in prevCont:
		tmpRow = [x for x in n]
		for m in prevCont:
			check = 0
			if m[2] != n[2] and m[1] == n[1]:
				tmpRow2 = [x for x in m]
				#n[len(n)-1] = float(n[len(n)-1]) + float(m[len(m)-1])
				sumValues = sumValues + float(m[len(m)-1])
				if subSpecies == []:
					tmpRow[0] = "Sub_Species"
					tmpRow[3] = ireplace(",", "", tmpRow[3])
					tmpRow[1] = ireplace("_complete_genome", "", tmpRow[3])
					subSpecies.append(tmpRow)
				else:
					tmpRow2[0] = "Sub_Species"
					tmpRow2[3] = ireplace(",", "", tmpRow2[3])
					tmpRow2[1] = ireplace("_complete_genome", "", tmpRow2[3])
					for t in subSpecies:
						if t[2] == tmpRow2[2]:
							check = 1
					if check != 1:
						subSpecies.append(tmpRow2)	
				check = 0
		if container == []:
			n[len(n)-1] = float(n[len(n)-1]) + sumValues
			container.append(n)
			sumValues = 0.0
		else:
			for c in container:
				if c[1] == n[1]:
					check = 1
					break
			if check != 1:
				n[len(n)-1] = float(n[len(n)-1]) + sumValues
				container.append(n)
				sumValues = 0.0
	for s in subSpecies:
		container.append(s)
	return container

def parseKraken2(fileInput):
	"""
	Columns in each row:
    	Column 1: percentage of reads in the clade/taxon in Column 6
    	Column 2: number of reads in the clade.
   		Column 3: number of reads in the clade but not further classified.
    	Column 4: code indicating the rank of the classification: 
									(U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, (S)pecies).
    	Column 5: NCBI taxonomy ID.
	
	Example run:
		['35.82', '1461', '1461', 'U', '0', 'unclassified']
		['64.18', '2618', '35', 'R', '1', 'root']
		['59.08', '2410', '11', 'R1', '131567', 'cellular', 'organisms']
		['58.69', '2394', '176', 'D', '2', 'Bacteria']
		['51.61', '2105', '85', 'P', '1224', 'Proteobacteria']
	
	Source: http://sepsis-omics.github.io/tutorials/modules/kraken/
	"""
	#Filter bacteria from results
	check_level = 1 # If level =1 stop appending to container
	
	addClasses = ['F','G','S']
	ClassLabel = ['Family', 'Genus', 'Species', 'Sub_Species', 'Sub_Genus', 'Sub_Family']
	#taxonLevel [Classification Char]
	taxonLevel = []
	#taxonName [Classification Level Name]
	taxonName = []
	
	container = []
	tmp = []
	tmpstr = ""
	run = 1
	f = open(fileInput)
	while run != 0:
		line = f.readline().split()
		if len(line) == 0:
			break
		"""
		FIX pattern some results have G1, S1,S2,S3 instead of G and S only
		regular expression match 
		pattern = re.compile("^([A-Z][0-9]+)+$")
		pattern.match(string)
		"""
		pattern = re.compile("D\Z")
		if pattern.match(str(line[3])):
			if "virus" not in str(line[5]) or "Virus" not in str(line[5]):
				check_level = 1
			if "virus" in str(line[5]) or "Virus" in str(line[5]):
				check_level = 0
		if check_level != 1:	
			if any(str(it) in str(line[3]) for it in addClasses):
				taxonLevel.append(line[3])
				tmp = line[5:len(line)]
				line = (line[0:5])
				for n in tmp:
					if tmpstr == "":
						tmpstr = tmpstr + n
					else:
						tmpstr = tmpstr + '_' + n 
				line.insert(0, tmpstr)
				taxonName.append(tmpstr)
				tmpstr = ""
				for n in range(len(addClasses)):
					pattern = re.compile("S[0-9]")
					if pattern.match(str(line[4])):
						line.insert(0, ClassLabel[3])
						container.append(line)
						break
					if str(line[4]) == str(addClasses[n]):
						line.insert(0, ClassLabel[n])
						container.append(line)
						break		   
				#print(line)
	###################################################
	tmpArr = []
	ListofArr = []
	tmpArrChar = []
	ListofArrChar = []
	for i in range(len(taxonLevel)):
		if i == len(taxonLevel)-1:
			ListofArr.append(tmpArr)
			ListofArrChar.append(tmpArrChar)
		if tmpArr == []:
			tmpArr.append(taxonName[i])
			tmpArrChar.append(taxonLevel[i])
		elif taxonLevel[i] == "F":
			ListofArr.append(tmpArr)
			tmpArr = []
			ListofArrChar.append(tmpArrChar)
			tmpArrChar = []
			tmpArr.append(taxonName[i])
			tmpArrChar.append(taxonLevel[i])
		else:
			tmpArr.append(taxonName[i])
			tmpArrChar.append(taxonLevel[i])		
	#########################################################
	tmpArr = []
	tmpArrChar = []
	check_level = 0
	for c in container:
		for i in range(len(ListofArrChar)):
			if c[1] in ListofArr[i] and check_level == 0:
				c.append(ListofArrChar[i])
				c.append(ListofArr[i])
				check_level = 1
		check_level = 0
	Genus = ""
	Species = ""
	######################################################
	for c in container:
		tmpArr = []
		check_level = 0
		for i in range(len(c[7])):
			if c[5] == "F":
				tmpArr.append(c[1])
				c[8] = taxonLevelString(tmpArr, 1)
				c.remove(c[7])
				break
			if c[5] == "G":
				if c[7][i] == "F":
					tmpArr.append(c[8][i])
					tmpArr.append(c[1])
					c[8] = taxonLevelString(tmpArr, 1)
					c.remove(c[7])
					break
			pattern = re.compile("S\Z")
			if pattern.match(c[5]):
				check_level = 1
				if c[7][i] == "F":
					tmpArr.append(c[8][i])
				if c[7][i] == "G":
					Genus = c[8][i]
				if c[8][i] == c[1]:
					tmpArr.append(Genus)
					tmpArr.append(c[1])
					c[8] = taxonLevelString(tmpArr, 1)
					c.remove(c[7])
					break
			pattern = re.compile("S[0-9]")
			if pattern.match(str(c[5])) and check_level == 0: 
				if c[7][i] == "F":
					tmpArr.append(c[8][i])
				if c[7][i] == "G":
					Genus = c[8][i]
				pattern = re.compile("S\Z")
				if pattern.match(c[7][i]):
					Species = c[8][i]
				if c[8][i] == c[1]:
					tmpArr.append(Genus)
					tmpArr.append(Species)
					tmpArr.append(c[1])
					c[8] = taxonLevelString(tmpArr, 1)
					c.remove(c[7])
					break
	return container

def combineResults(Kraken, FVE, Meta, ID):
	"""
	#Create header 
	# Patient_ID	Name	Kraken	FastViromeExplorer	MetaPhlan	Classification	Full_Taxon_level
	Header = "Patient_ID\tName\tKraken_#nrReads\tFastViromeExplorer\tMetaPhlan\tClassification\tFull_Taxon_level(Family->Genus->Species->Sub_species)\n"
	CreateFile(ID)
	name = "Combined_Results_"+ID+".txt"
	f=open(name, "a+")
	f.write(Header)
	f.close() 
	"""
	#tmp = [ID,"Name","Kraken", "FVE", "Meta", "Class", "level"]
	tmp = [ID,"-","-", "-", "-", "-", "-"]
	container = []
	check = 1
	for k in Kraken:
		#Name
		tmp[1] = k[1]
		#Class
		tmp[5] = k[0]
		#Kraken
		tmp[2] = k[3]
		#Taxon level from Kraken
		tmp[6] = k[7]
		for f in FVE:
			if tmp[1] == f[1] and tmp[5] == f[0]:
				tmp[3] = f[len(f)-1]
				#tmp[6] = f[4].replace("Unclassified;", "")
				#tmp[6] = tmp[6].replace(";","->")
			for m in Meta:
				if tmp[1] == m[1] and tmp[5] == m[0]:
					tmp[4] = m[len(m)-1]
		container.append(tmp)
		tmp = [ID,"-","-", "-", "-", "-", "-"]
		
	for f in FVE:
		tmp = [ID,"-","-", "-", "-", "-", "-"]
		check = 0
		for k in Kraken:
			if f[1] == k[1] and f[0] == k[0]:
				check = 1
		if check != 1:
			#Name
			tmp[1] = f[1]
			#Class
			tmp[5] = f[0]
			#FVE
			tmp[3] = f[len(f)-1]
			#level
			tmp[6] = f[4].replace("Unclassified;", "")
			tmp[6] = tmp[6].replace(";","->")
			for m in Meta:
				if tmp[1] == m[1] and tmp[5] == m[0]:
					tmp[4] = m[len(m)-1]
			container.append(tmp)
		
	for m in Meta:
		tmp = [ID,"-","-", "-", "-", "-", "-"]
		check = 0
		for f in FVE:
			if f[1] == m[1] and f[0] == m[0]:
				check = 1
			for k in Kraken:
				if k[1] == m[1] and k[0] == m[0]:
					check = 1
		if check != 1:
			#Name
			tmp[1] = m[1]
			#Class
			tmp[5] = m[0]
			#Meta
			tmp[4] = m[len(m)-1]
			tmp[6] = taxonLevelString(m[2].split("|"), 2)
			container.append(tmp)
	#for c in container:
	#	print(c)
	WriteToFile(container, ID)
	return None

def main():
    #newFileName = "LCH_episome_mapped_reads.tsv"
    #files = glob.glob('/proj/virus/results/LCH/190405/bbmap_results/*/*_flagstat.txt')
    #newFileName = sys.argv[2]
    #files = glob.glob(str(sys.argv[1]+"*_flagstat.txt"))
    #files = glob.glob(str(sys.argv[1]))
	
	#PATHWAY
	files_metaphlan2 = glob.glob('/home/amanj/Documents/ttv_bonk/20190401/discovery/P*/reads/*metaphlan2.tsv')
	files_fastviromeexplorer = glob.glob('/home/amanj/Documents/ttv_bonk/20190401/discovery/P*/reads/*fastviromeexplorer_abundance.tsv')
	files_kraken2 = glob.glob('/home/amanj/Documents/ttv_bonk/20190401/discovery/P*/reads/*kraken2_report.txt')
	files_metaphlan2.sort()
	files_fastviromeexplorer.sort()
	files_kraken2.sort()
	
	for i in range(len(files_metaphlan2)):
		fileName = files_metaphlan2[i]
		M_container = parseMetaphlan2(fileName)
		#Fetching Patient ID
		Patient_ID = M_container[0][0]
		M_container.pop(0)
		fileName = files_fastviromeexplorer[i]
		#fileName = "/home/amanj/Documents/190425/pcr_products_180912/discovery/P11463_1001_S1_L001/reads/P11463_1001_S1_L001_fastviromeexplorer_abundance.tsv"
		F_container = parseFVE(fileName)
		F_container = FVEaddSubSpecies(F_container)
		fileName = files_kraken2[i]
		K_container = parseKraken2(fileName)
		combineResults(K_container, F_container, M_container, Patient_ID)
	
	all_results = glob.glob('Combined_Results_*.txt')
	all_results.sort()
	
	#Create header 
	# Patient_ID	Name	Kraken	FastViromeExplorer	MetaPhlan	Classification	Full_Taxon_level
	Header = "Patient_ID\tName\tKraken_nrReads\tFastViromeExplorer\tMetaPhlan\tClassification\tFull_Taxon_level(Family-Genus-Species-Sub_species)\n"
	f = open("All_patient_results.txt", "w")
	f.write(Header)
	for result in all_results:
		f1 = open(result, "r")
		f.write(f1.read())
		f1.close()
	f.close()
	"""
	f1 = open("All_patient_results_into_one_file.txt", "r")
	titles = [string.strip(t) for t in string.split(f1.readline(), sep="\t")]
	for l in sys.stdin:
		d = {}
		for t, f in zip(titles, string.split(l, sep="\t")):
			d[t] = f
		print(common.json.dumps(d, indent=4))
	"""
	return None




main()