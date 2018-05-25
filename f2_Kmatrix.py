#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# f2_Kmatrix.py
# Created: 2018-05-08 17:59:31.697511
# Author: Fabian Krinner
import os, sys
from KmatrixFitResults import KmatrixFitResult
import ROOT
from time import sleep
import scipy.optimize
from rootfabi import root_open
from globalDefinitions import mPi
def fPlot(v):
	return v.imag

def findPole(Kmatrix, startValues, maxPolVal = 1.e-8):

	res = scipy.optimize.minimize(Kmatrix.absInverse,startValues)
	if not abs(res.fun) < maxPolVal:
		return None
	return res.x

def findAllPoles(Kmatrix, reRange, imRange, nSteps = 4, maxDist = 1.e-5):
	stepRe = (reRange[1] - reRange[0])/nSteps
	stepIm = (imRange[1] - imRange[0])/nSteps
	poles  = []
	for i in range(nSteps):
		re = reRange[0] + (i+.5)*stepRe
		for j in range(nSteps):
			im   = imRange[0] + (j+.5)*stepIm
			pole = findPole(Kmatrix, [re,im])
			if pole is None:
				continue
			old = False
			for p in poles:
				dist = abs(p[0] - pole[0])**2 + abs(p[1] - pole[1])**2
				if dist < maxDist:
					old = True
					break
			
			if not old:
				poles.append(pole)
	return poles

nBinsPlot  = 250
reMin      =-  .25
reMax      =  6.
imMin      =- 2.
imMax      =  2.


histAllPolePos     = ROOT.TH2D("hHh","hHh",         nBinsPlot, reMin, reMax ,nBinsPlot, 0.,imMax)
histAllPolePos_CM  = ROOT.TH2D("hHh_CM","hHh_CM",   nBinsPlot, reMin, reMax ,nBinsPlot, 0.,imMax)
histAllPolePos_rho = ROOT.TH2D("hHh_rho","hHh_rho", nBinsPlot, reMin, reMax ,nBinsPlot, 0.,imMax)

def domtBin(tBin, mBin):
	ROOT.gStyle.SetOptStat(0)
	findBest     = True
	allFileNames = []
	if findBest:
		tString   = "_t"+str(tBin)+"_"
		bestChi2  = float('inf')
		bestFn    = None
		worstChi2 =-float('inf')
#		folder    = "./KmatrixResults/"
		folder    = "./KmatrixRestrictedZMfixing_allFourSectorsActive/"
		m3Pi      = .5 + (mBin+.5)*.04
		for fn in os.listdir(folder):
			if not tString in fn:
				continue
			if not "_kPol0_" in fn:
				continue
			if not "m"+str(mBin)+"-"+str(mBin+1) in fn:
				continue
			if not fn.startswith("2mpF2_"):
				continue
#			if not "_nPol3_" in fn:
#				continue
#			if not "_CM_" in fn:
#				continue
			fileName     = folder + fn
			result       = KmatrixFitResult(fileName)
			try:
				result.chi2
				allFileNames.append((result.chi2, fileName,result.chi2/result.ndf))
			except AttributeError:
				print "ERROR in: ",fn
				continue
			worstChi2 = max(worstChi2, result.chi2)
			if result.chi2 < bestChi2:
				bestChi2 = result.chi2 
				bestFn   = fn
			try:
				inFileName = folder + bestFn
			except:
				print "ERROR : : : : : : "
				continue
	else:
		inFileName = "./KmatrixResults/2mpF2_Kmatrix_nPol3_rho_kPol2_pPol7-4_t0_m40-50_9696.dat"


	allFileNames.sort()
	allFileNames.reverse()
	allFileNames = allFileNames[~1:]



	c2h_rho = ROOT.TH1D("c2","c2",250, bestChi2, worstChi2)
	c2h_CM  = ROOT.TH1D("c2","c2",250, bestChi2, worstChi2)
	c2h_CM.SetLineColor(2)


	first = True
	hsh = []
	for f,inFileName in enumerate(allFileNames):
#		if not f == 182:
#			continue
		print f+1,'/',len(allFileNames), inFileName
		result     = KmatrixFitResult(inFileName[1])
		pyCommand = "python 2mp_data_zeroModeFixingComparison.py "+str(result.tBin)+ " ifn="+result.inFileName+ " mBins-"+str(result.mRange[0])+'-'+str(result.mRange[1])+" "+result.phaseSpace+" m3pol"+str(result.pPol[0])+" m2pol"+str(result.pPol[1])+" pdpo"+str(result.kPol+1)
		print ">>>",pyCommand,"<<<"
		Kmatrix    = result.produceFunction()
		pVector    = result.producePvector()

		intensHist = ROOT.TH1D("ntns",str(tBin)+"_"+str(mBin), nBinsPlot, 2*mPi, m3Pi - mPi)

		masses = [intensHist.GetXaxis().GetBinCenter(i+1) for i in range(nBinsPlot)]
		vK = Kmatrix(masses)
		vP = pVector(masses,[m3Pi])
		for i in range(nBinsPlot):
			intensHist.SetBinContent(i+1, abs(vK[i]*vP[i])**2)
#			intensHist.SetBinContent(i+1, vP[i])

		mMinNorm   = (2*mPi)
		mMaxNorm   = (m3Pi-mPi)
		binMinNorm = intensHist.GetXaxis().FindBin(mMinNorm)
		binMaxNorm = intensHist.GetXaxis().FindBin(mMaxNorm)
		norm       = 0.
		for i in range(binMinNorm, binMaxNorm):
			norm += intensHist.GetBinContent(i)
		intensHist.Scale(1./norm)

		hsh.append(intensHist)

		hist = ROOT.TH2D("hhh",str(tBin)+"_"+str(mBin), nBinsPlot, reMin, reMax ,nBinsPlot, imMin,imMax)

		BWtoSet = [(1.272,0.1867),(1.525, 0.075),(2.011,0.202)]

		sToSet = [BW[0]**2 for BW in BWtoSet]
		GtoSet = [BW[0]*BW[1] for BW in BWtoSet]

		bwSetLength   = 0.1
		setVal        = 10.
	
#		for iX in range(nBinsPlot):
#			x = hist.GetXaxis().GetBinCenter(iX+1)
#			for iY in range(nBinsPlot):
#				y = hist.GetYaxis().GetBinCenter(iY+1)			
#				s = x+1.j*y
#				val = Kmatrix.complexCall(s)
#				hist.SetBinContent(iX+1, iY+1,fPlot(val))
#
#		for s in sToSet:
#			iX = hist.GetXaxis().FindBin(s)
#			for iY in range(nBinsPlot):
#				hist.SetBinContent(iX, iY+1, setVal)
#	
#		hist.SetMaximum( 10.)
#		hist.SetMinimum(-10.)
#	
#		hist.Draw("COLZ")
#		raw_input("press <enter> to go to the second sheet")
		Kmatrix.secondSheet = True
		allPoles = findAllPoles(Kmatrix, [-1.,6.], [0.,2.])
		for p in allPoles:
			histAllPolePos.Fill(p[0],abs(p[1]))
			if "_CM_" in inFileName[1]:
				histAllPolePos_CM.Fill(p[0],abs(p[1]))
			else:
				histAllPolePos_rho.Fill(p[0],abs(p[1]))
		for iX in range(nBinsPlot):
			x = hist.GetXaxis().GetBinCenter(iX+1)
			for iY in range(nBinsPlot):
				y = hist.GetYaxis().GetBinCenter(iY+1)			
				s = x+1.j*y
				val = Kmatrix.complexCall(s)
				hist.SetBinContent(iX+1, iY+1, fPlot(val))
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
		hist.SetMaximum( 10.)
		hist.SetMinimum(-10.)
		if not "_CM_" in inFileName[1]:
			c2h_rho.Fill(inFileName[0])
		else:
			c2h_CM.Fill(inFileName[0])

		c2.cd()
#		c2h_rho.Draw()
#		c2h_CM.Draw()
#		histAllPolePos.Draw("COLZ")
		if first:
			intensHist.Draw()
		else:
			intensHist.Draw("SAME")
		c2.Update()
		c1.cd()
		hist.Draw("COLZ")
		c1.Update()
#		raw_input("press <enter> to quit")
#		sleep(.3)

		first = False
#	with root_open("./poleDistribution.root","RECREATE"):
#		c2h_rho.Write()
#		c2h_CM.Write()
#		histAllPolePos.Write()
#		histAllPolePos_CM.Write()
#		histAllPolePos_rho.Write()
#	raw_input("<enter> to exit")
#	c1.Close()
#	c2.Close()
#	c1.Delete()
#	c2.Delete()
#	raw_input()

c1      = ROOT.TCanvas()
c2      = ROOT.TCanvas()

def main():

	for t in range(4):
		for b in range(25, 40):
			domtBin(t,b)
	c1.cd()
	histAllPolePos.Draw("COLZ")
	c1.Update()
	raw_input()

if __name__ == "__main__":
	main()
