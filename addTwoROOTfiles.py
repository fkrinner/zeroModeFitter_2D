# addTwoROOTfiles.py
import os, sys
from rootfabi import root_open, GetKeyNames

def getHistograms(inFileName):
	hists = {}
	with root_open(inFileName, "READ") as inFile:
		names = GetKeyNames(inFile)
		for nn in names:
			hist = inFile.Get(nn)
			hist.SetDirectory(0)
			hists[nn] = hist
	return hists

def main():
	inFile1 = "980.root"
	inFile2 = "1500.root"
	outFile = "f0s.root"
	hists1  = getHistograms(inFile1)
	hists2  = getHistograms(inFile2)
	for nn in hists1:
		hists1[nn].Add(hists2[nn])
	with root_open(outFile, "RECREATE"):
		for nn in hists1:
			hists1[nn].Write()

if __name__ == "__main__":
	main()
