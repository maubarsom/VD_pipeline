import sys
import os

def WriteToFile(nr_mapped, proc_mapped, FileType, FileSampleID, Ref, name):
	if os.path.exists(sys.argv[2]):

	else:
		str1 = "SampleID" + "\t" + "Reference" + "\t" + "Type" + "\t"  + "n_mapped"  + "\t" + "pct_mapped" + "\n"
		f=open(name, "a+")
		f.write(str1)
		f.close()
		print("else")

	str1 = FileSampleID + "\t" + Ref + "\t" + FileType + "\t" + nr_mapped + "\t" + proc_mapped + "\n"
	f=open(name, "a+")
	f.write(str1)
	f.close()
	return None

def main():
	tmpLine = []
	f = open(sys.argv[1])
	line = f.readline()
	tmpLine = line.split()
	newFileName = sys.argv[2]
	FileName = sys.argv[1]
	FileName = FileName.split("_")
	fileType = FileName[3]
	SampleID = FileName[0] + "_" + FileName[1]
	ref = FileName[2]
	while line:
		line = line.strip()
		line = f.readline()
		tmpLine = line.split()
		if tmpLine:
			if tmpLine[3] == "mapped":
				WriteToFile(tmpLine[0], tmpLine[4] + tmpLine[5] + tmpLine[6], fileType, SampleID, ref, newFileName)
		else:
			break
	f.close()
	return None
main()

