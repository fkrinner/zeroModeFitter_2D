#!/usr/bin/python
# analyze_vff.py
# Created: 2018-07-26 17:55:07.998683
# Author: Fabian Krinner
import os, sys
import re
import ROOT

def makeValHists(valMap, nBins = 100, delta = 1.e-5):
	histMap = {}
	for v in valMap:
		valMap[v].sort()
		if v == "c2ndf":
			hist = ROOT.TH1D(v,v,nBins, 0., valMap[v][~0]+delta)
		else:
			hist = ROOT.TH1D(v,v,nBins, valMap[v][0], valMap[v][~0]+delta)
		for p in valMap[v]:
			hist.Fill(p)
		histMap[v] = hist
	return histMap

class vff_result:
	def __init__(self, inFileName):
		self.inFileName = inFileName
		self.numLim = 1.e-10
		self.readFile()

	def readFile(self):
		self.params = {}
		with open(self.inFileName, 'r') as inFile:
			inPar = False
			for line in inFile.readlines():
				if "- - - - parameters - - - -" in line:
					inPar = True
					continue
				if "- - - - - fit - - - - -" in line:
					inPar = False
					continue		
				if inPar:
					chunks = line.split()
					self.params[chunks[0]] = (float(chunks[2]),float(chunks[4]))
					continue
				if "chi2/NDF:" in line:
					chunks = re.split(' |/|=|\n',line)
					self.chi2  = float(chunks[2])
					self.ndf   = int(chunks[3])
					self.c2ndf = float(chunks[4])
					if abs(self.chi2/self.ndf - self.c2ndf) > self.numLim:
						print self.chi2, self.ndf, self.chi2/self.ndf, self.c2ndf,self.chi2/self.ndf-self.c2ndf
						raise ValueError("c2ndf mismatch")

def main():
	folderName = "./vectorFormFactorResults_1mp/"

	valMap = {'c2ndf':[]}

	for fn in os.listdir(folderName):
		res = vff_result(folderName + fn)
		valMap['c2ndf'].append(res.c2ndf)
		for v in res.params:
			if not v in valMap:
				valMap[v] = []
			valMap[v].append(res.params[v][0])
	histMap = makeValHists(valMap)
	rightParAlignment = ['c2ndf','lambP','lambPP','M','G','MP','GP','MPP','GPP','alP','phP','alPP','phPP']
	c1 = ROOT.TCanvas()
	for v in rightParAlignment:
		histMap[v].Draw()
		c1.Update()
		raw_input()


if __name__ == "__main__":
	main()
