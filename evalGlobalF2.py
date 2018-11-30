#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# evalGlobalF2.py
# Created: 2018-05-28 11:01:23.898231
# Author: Fabian Krinner
import os, sys
from parameterTrackingParameterizations import parameter, simpleOneChannelKmatrix
from globalDefinitions import mPi
import ROOT

from rootfabi import root_open

def makeKmatrix(inFileName):
	nPol  = None
	kPol  = None
	useCM = None
	for a in inFileName.split("_"):
		if a.startswith("nPol"):
			nPol = int(a[4:])
		if a.startswith("kPol"):
			kPol = int(a[4:])
		if a == "CM":
			useCM = True
		if a == "rho":
			useCM = False

	if nPol is None:
		raise IOError("Could not determine nPol")
	if kPol is None:
		raise IOError("Could not determine kPol")
	if useCM is None:
		raise IOError("Could not determine useCM")
	parameters = []
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			if "- - - " in line:
				continue
			if "t" in line:
				break
			chunks    = line.split()
			par       = parameter(float(chunks[2]),chunks[0], float(chunks[4]))
			parameters.append(par)
	if not len(parameters) == 2*nPol + kPol+1:
		print nPol,kPol+1,"=>",2*nPol + kPol+1,"!=",len(parameters)
		raise ValueError("Number of parameters does not match")
	retVal = simpleOneChannelKmatrix(parameters, nPol, kPol+1,4*mPi**2)
	retVal.use_CM = useCM

	return retVal

def readChi2ndf(inFileName):
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			if "chi2/NDF:" in line:
				chi2 = float(line.split('=')[-1])
				return chi2

def main():
	ROOT.gStyle.SetOptStat(0)
	folder   = "./globalKmatrixFits/"
	bestChi2 = float("inf")
	bestFile = None
	c2s = []
	for fn in os.listdir(folder):
#		if not "_polyDegP3_" in fn:
#			continue
		chi2 = readChi2ndf(folder + fn)
		c2s.append(chi2)
		print " > > > >",chi2,"< < < < <"
		if chi2 < bestChi2:
			bestChi2 = chi2
			bestFile = fn
	print "Best chi2 in",bestChi2, "from", bestFile
	c2s.sort()
	histC2 = ROOT.TH1D("chi2", "chi2", 1000, c2s[0], c2s[~0])
	for v in c2s:
		histC2.Fill(v)
	histC2.Draw()
	raw_input()

	inFileName = folder + bestFile
	Kmatrix    = makeKmatrix(inFileName)
	Kmatrix.secondSheet = True
	nBinsPlot  = 500
	reMin      =-  .25
	reMax      =  6.
	imMin      =- 2.
	imMax      =  2.
	BWtoSet = [(1.272,0.1867),(1.525, 0.075),(2.011,0.202)]

	BWtoSet = []

	sToSet = [BW[0]**2 for BW in BWtoSet]
	GtoSet = [BW[0]*BW[1] for BW in BWtoSet]

	bwSetLength   = 0.1
	setVal        = 10.



	hist = ROOT.TH2D("hhh","Kmatrix_second_sheet", nBinsPlot, reMin, reMax ,nBinsPlot, imMin,imMax)
	for iX in range(nBinsPlot):
		x = hist.GetXaxis().GetBinCenter(iX+1)
		for iY in range(nBinsPlot):
			y = hist.GetYaxis().GetBinCenter(iY+1)			
			s = x+1.j*y
			val = Kmatrix.complexCall(s)
			hist.SetBinContent(iX+1, iY+1, val.imag)
	for s in sToSet:
		iX = hist.GetXaxis().FindBin(s)
		for iY in range(nBinsPlot):
			hist.SetBinContent(iX, iY+1, -setVal)
		for g,im in enumerate(GtoSet):
			re = sToSet[g]
			binMin = hist.GetXaxis().FindBin(re-bwSetLength)
			binMax = hist.GetXaxis().FindBin(re+bwSetLength)
			binIm  = hist.GetYaxis().FindBin(im)
			binImM = hist.GetYaxis().FindBin(-im)
			for b in range(binMin, binMax):
					hist.SetBinContent(b+1, binIm, -setVal)
					hist.SetBinContent(b+1, binImM, -setVal)

	with root_open("Kmatrix.root","RECREATE"):
		hist.Write()

	hist.Draw("colz")
	raw_input()

if __name__ == "__main__":
	main()
