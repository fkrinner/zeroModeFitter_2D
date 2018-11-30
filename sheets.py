#!/usr/bin/python
# sheets.py
# Created: 2018-08-20 13:20:27.004587
# Author: Fabian Krinner
import os, sys
from physUtils import chewMandelstam as chms
from cmath import exp, log, pi

import scipy.optimize as opt

import ROOT

class thresholdFunction:
	def setSheet(self, sheet = 0):
		if sheet > self.maxSheet:
			raise ValueError(self.name + " has a maxSheet of " + str(self.maxSheet))
		self.sheet = sheet		

class iSigma(thresholdFunction):
	def __init__(self, sThresh, name):
		self.name     = name
		self.sThresh  = sThresh
		self.sheet    = 0
		self.maxSheet = 1

	def __call__(self, s):
		retVal = -(self.sThresh/s-1.+0.j)**.5
		if self.sheet == 0:
			return -retVal
		else:
			return retVal

class chewMandelstam(thresholdFunction):
	def __init__(self, m1, m2, name):
		self.name     = name
		self.m1       = m1
		self.m2       = m2
		self.sThresh  = (m1+m2)**2
		self.sheet    = 0
		self.maxSheet = 1
	
	def __call__(self, s):
		retVal = -chms(s,self.m1, self.m2)
		if self.sheet == 0:
			return retVal
		else:
			return retVal - 2*(self.sThresh/s-1.+0.j)**.5

class PwavePhaseSpace(thresholdFunction):
	def __init__(self, mP, mu2, name):
		self.name     = name
		self.mP       = mP
		self.sThresh  = 4*mP**2
		self.sig      = iSigma(self.sThresh, name+"_sigma")
		self.mu2      = mu2
		self.sheet    = 0
		self.maxSheet = 1

	def __call__(self,s):
		sig    = self.sig(s)/1.j
		retVal = log(self.mP**2/self.mu2) + 8*self.mP**2/s - 5./3. + sig**3 * log((sig+1.)/(sig-1.))
		if self.sheet == 0:
			return retVal
		else:
			return retVal + 2*pi* (self.sThresh/s-1.+0.j)**1.5

def doCircularPlot(funcs, absS, nBins = 1000):
	if len(funcs) == 0:
		raise RuntimeError("doCircularPlot: No functions given")
	histsReal = []
	histsImag = []
	histZL    = ROOT.TH1D("zeroLine","zeroLine",nBins,-pi,pi)
	histZL.SetLineStyle(2)
	for f in range(len(funcs)):
		histsReal.append(ROOT.TH1D("real_"+str(f),"real_"+str(f),nBins,-pi,pi))
		histsReal[f].SetLineStyle(9)
		histsImag.append(ROOT.TH1D("imag_"+str(f),"imag_"+str(f),nBins,-pi,pi))
	for b in range(nBins):
		phase = histsReal[0].GetXaxis().GetBinCenter(b+1)
		s = absS*exp(1.j*phase)
		for f, func in enumerate(funcs):
			val = func(s)
			histsReal[f].SetBinContent(b+1, val.real)
			histsImag[f].SetBinContent(b+1, val.imag)
	for f in range(1,len(funcs)):
		histsReal[f].SetLineColor(f+1)
		histsImag[f].SetLineColor(f+1)
	maxx = -float('inf')
	minn =  float('inf')
	for f in range(len(funcs)):
		maxx = max(maxx, histsReal[f].GetMaximum(), histsImag[f].GetMaximum())
		minn = min(minn, histsReal[f].GetMinimum(), histsImag[f].GetMinimum())
	histsImag[0].SetMaximum(maxx)
	histsImag[0].SetMinimum(minn)

	histsImag[0].Draw()
	histsReal[0].Draw("SAME")
	for f in range(1, len(funcs)):
		histsImag[f].Draw("SAME")
		histsReal[f].Draw("SAME")
	histZL.Draw("SAME")
	raw_input()

def findPoles(func, reMin, reMax, imMin, imMax, nSteps = 10, delta = 1.e-4):
	def minF(p):
		s = p[0] + 1.j*p[1]
		return 1./abs(func(s))**2
	reStep = (reMax-reMin)/nSteps
	imStep = (imMax-imMin)/nSteps
	poles = []
	for i in range(nSteps):
		reStart = reMin + (i+.5)*reStep
		for j in range(nSteps):
			imStart = imMin + (j+.5)*imStep
			res = opt.minimize(minF, [reStart, imStart])
			if res.fun < delta:
				poles.append(res.x[0] + 1.j * res.x[1])
	poles = checkAndIdentifyPoles(func, poles, delta)
	return poles

def checkAndIdentifyPoles(func, poles, delta, onlyUpperHalfPlane = True):
	retVal = []
	for p in poles:
		if p.imag < 0. and onlyUpperHalfPlane:
			continue
		plusRe = func(p + delta    )
		minuRe = func(p - delta    )
		plusIm = func(p + 1.j*delta)
		minuIm = func(p - 1.j*delta)
#		print plusRe/minuRe
#		print plusIm/minuIm

#		if not plusRe.real*minuRe.real * plusRe.imag*minuRe.imag < 0.:
#			print plusRe, minuRe
#			raise ValueError("No sing change for pole "+str(p)+" in real direction")
#		if not plusIm.real*minuIm.real* plusIm.imag*minuIm.imag < 0.:
#			print plusIm,minusIm
#			raise ValueError("No sing change for pole "+str(p)+" in imag direction")
		found = False
		for q in retVal:
			if abs(p - q) < delta:
				found = True 
				break
		if not found:
			retVal.append(p)	
	return retVal

def main():
	m1 = .1
	m2 = .1

	sTh = (m1+m2)**2

	cm    = chewMandelstam(m1, m2, "tschu_mandelschtamm")
	cmss  = chewMandelstam(m1, m2, "tschu_mandelschtamm_secondSheet")
	cmss.setSheet(1)

	sig   = iSigma((m1+m2)**2, "sigmer")
	sigss = iSigma((m1+m2)**2, "sigmer_secondSheet")
	sigss.setSheet(1)

	Pps   = PwavePhaseSpace((m1+m2)/2, 1.,"P_face_spase")

	Ppsss = PwavePhaseSpace((m1+m2)/2, 1.,"P_face_spase_secondSheet")
	Ppsss.setSheet(1)

	doCircularPlot([cm, cmss ], sTh+0.01)
	doCircularPlot([Pps,Ppsss], sTh+0.01)


if __name__ == "__main__":
	main()
