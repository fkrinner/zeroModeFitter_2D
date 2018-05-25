# compareTheo.py
# Created: 2018-03-06 13:59:35.267712
# Author: Fabian Krinner
import os, sys

def loadTxtFile(inFileName, reIndex = 0, d = 1):
	data = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			data.append(float(chunks[reIndex]) + 1.j*float(chunks[reIndex+d]))
	return data

def main():
	inFileMasses = "./theo_3pp0p2pp_32.intens"
	masses       = loadTxtFile(inFileMasses, 0)

#	inFileFabi = "./3pp0p2pp_theoAmplitude_deleteMeSoon.dat"
	inFileFabi = "./theo_3pp0p2pp_32.argand"
	dataFabi   = loadTxtFile(inFileFabi, 0, 2)

	inFileDima = "/nfs/freenas/tuph/e18/project/compass/analysis/dryabchi/free_isob_3pp/0.100-0.141/param_14_01.dat"
	dataDima   = loadTxtFile(inFileDima, 7, 1)

	sumf = 0.
	sumd = 0.
	for i in range(min(len(dataFabi), len(dataDima))):
		sumf += abs(dataFabi[i])
		sumd += abs(dataDima[i])
	

	with open("./ratioFile_deleteMe.dat", 'w') as outFile:
		for i in range(min(len(dataFabi), len(dataDima))):
			outFile.write(str(masses[i].real)+" "+str(abs(dataFabi[i])/(abs(dataDima[i])*sumf/sumd))+" "+str((abs(dataDima[i])*sumf/sumd))+"\n")

if __name__ == "__main__":
	main()
