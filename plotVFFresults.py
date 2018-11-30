#!/usr/bin/python
# plotVFFresults.py
# Created: 2018-08-21 17:49:53.600052
# Author: Fabian Krinner
import os, sys
import ROOT

ndfs = [52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 96, 96]
def parseFile(inFileName):
	vals = []
	nBin = int(inFileName.split("_")[~1][:2])
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			if len(chunks) == 2:
				vals.append(float(chunks[1]))
			if len(chunks) == 1:
				cchunks = chunks[0].split('/')
				chi2 = float(cchunks[0])
				if len(cchunks) == 0:
					ndf  = ndfs[nBin-25]
					print "Using default NDF of",ndf
				else:
					ndf = int(cchunks[1])
				vals.append(chi2/ndf)
	return vals

def getNames(inFileName):
	names = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			if len(chunks) == 2:
				names.append(chunks[0])
	names.append("chi2/ndf")	
	return names

def main():
	inFolderName = "./global_vff_fits_subThresh/"
	values       = []
	names        = None
	for fn in os.listdir(inFolderName):
		if names is None:
			names = getNames(inFolderName + fn)
		values.append(parseFile(inFolderName + fn))
	
	nBins = 20
	c1    = ROOT.TCanvas()

	for i in range(len(values[0])):
		minn = min([v[i] for v in values])
		maxx = max([v[i] for v in values])
		if i == len(values[0])-1:
			minn = 0.
			maxx = 20.

		hist = ROOT.TH1D(names[i], names[i], nBins, minn, maxx)
		for v in values:
			hist.Fill(v[i])
		hist.Draw()
		c1.Update()
		raw_input(names[i])

if __name__ == "__main__":
	main()
