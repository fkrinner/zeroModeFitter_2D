#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
# KmatrixFitResults.py
# Created: 2018-04-13 16:30:27.867556
# Author: Fabian Krinner
import os, sys
from parameterTrackingParameterizations import parameter, simpleOneChannelKmatrix, twoDimensionalRealPolynomial
from globalDefinitions import mPi

from math import sin,cos,pi

import ROOT

knownWaves = ['1pp','1mp','2mpF2']
class KmatrixFitResult:
	def __init__(self, inFileName, numLim = 1.e-13):
		self.inFileName = inFileName
		self.numLim     = numLim
		fileName        = inFileName.split(os.sep)[-1]
		chunks          = fileName.split('_')
		self.wave       = None
		self.mRange     = None
		self.tBin       = None
		self.phaseSpace = None
		self.nPol       = None
		self.kPol       = None
		self.pPol       = None
		self.range      = None
		self.coma       = None
		self.BWresult   = None
		self.sThresh    = 4*mPi**2
		self.nCmx       = 0
		for c in chunks:
			if c.endswith(".dat"):
				c = c[:-4]
			if c == 'Kmatrix':
				pass
			elif c in knownWaves:
				if self.wave is None:
					self.wave = c
				else:
					raise IOError("Wave set twice")
			elif c in ['rho','CM']:
				if self.phaseSpace is None:
					self.phaseSpace = c
				else:
					raise IOError("phase space set twice")
			elif c.startswith('range'):
				if self.range is None:
					self.range = float(c[5:])
				else:
					raise IOError("range set twice")
			elif c.startswith('nPol'):
				if self.nPol is None:
					self.nPol = int(c[4:])
				else:
					raise IOError("nPol set twice")
			elif c.startswith("nCmx"):
				self.nCmx = int(c[4:])

			elif c.startswith('kPol'):
				if self.kPol is None:
					self.kPol = int(c[4:])
				else:
					raise IOError("kPol set twice")
			elif c.startswith('pPol'):
				if self.pPol is None:
					self.pPol = (int(c[4:].split('-')[0]),int(c[4:].split('-')[1]))
				else:
					raise IOError("kPol set twice")
			elif c.startswith('t'):
				if self.tBin is None:
					self.tBin = int(c[1:])
				else:
					raise IOError("tBin set twice")
			elif c.startswith('m'):
				if self.mRange is None:
					self.mRange = (int(c[1:].split('-')[0]), int(c[1:].split('-')[1]))
				else:
					raise IOError("mRange set twice")
			else:
				try:
					int(c)
				except ValueError:
					raise IOError("Unparsable component of file name: '" + c + "'")
		self.loadFile()

	def loadFile(self):
		parameters = []
		with open(self.inFileName) as inFile:
			inPars = False
			for line in inFile.readlines():
				if "parameters" in line:
					inPars = True
					continue
				if "fit" in line and not "_fit_" in line:
					inPars = False
					continue
				if inPars:
					chunks    = line.split()
					par       = parameter(float(chunks[2]),chunks[0], float(chunks[4]))
					parameters.append(par)
				if "chi2/NDF" in line:
					self.chi2 = complex(line.split()[1].split('/')[0]).real
					self.ndf  = int(line.split('/')[2].split('=')[0])
				if "BW par" in line:
					chunks = line.split()
					self.mRho = float(chunks[2])
					self.Grho = float(chunks[3])
					if len(chunks) > 7:
						fVal      = float(chunks[7])
						if abs(fVal) > self.numLim:
							raise ValueError("Rho pole finding did not converge!!! Function value: "+str(fVal))
						self.fVal = fVal
				if "coma" in line:
					if self.coma is None:
						exec line.replace("coma","self.coma = ")	

					else:
						raise ValueError("COMA set twoce")

				if "BW' par" in line:
					chunks         = line.split()
					self.mRhoPrime = float(chunks[2])
					self.GrhoPrime = float(chunks[3])
					if len(chunks) > 7:
						self.fValPrime = float(chunks[7])
				if "BW_fit_result" in line:
					if self.BWresult is None:
						chunks = line.split()
						vals = [float(chunks[i]) for i in range(1, 5)]
						self.BWresult = ((vals[0],vals[1]),(vals[2],vals[3]))
					else:
						raise ValueError("BWpar set twice")
		self.parameters = parameters
		self.nPar       = len(parameters)
		if self.nPol is None:
			self.nPol       = (self.nPar - (self.kPol+1) - (self.pPol[0]+1) *self.pPol[1])/2	
			if not self.nPol == 1:
				print self.inFileName
				raise NotImplementedError("Number of poles not implemented yet")
	#	if not self.checkPole():
	#		raise ValueError("Invalid pole position encountered")

	def makeMassWidthComa(self):
		if self.coma is None:
			raise RuntimeError("No coma present")
		jacobian = [[1./(2*self.mRho),0.],[-self.Grho,-1./self.mRho]]
		mGcoma   = [[0.,0.],[0.,0.]]
		for i in range(2):
			for j in range(2):
				for k in range(2):
					for l in range(2):
						mGcoma[i][l] += self.coma[j][k]*jacobian[i][j]*jacobian[l][k]
		return mGcoma

	def producePvector(self):
		nPar    = 2*self.nPol + self.kPol + 1 + 3*self.nCmx
		pVector = twoDimensionalRealPolynomial(self.pPol[1],self.pPol[0], self.parameters[nPar:], baseExponent = 2)
		return pVector

	def produceFunction(self):
		nPar   = 2*self.nPol  + self.kPol +  1 + 3*self.nCmx
		retVal = simpleOneChannelKmatrix(self.parameters[:nPar], self.nPol, self.kPol+1, self.sThresh, nComplex = self.nCmx)
		retVal.use_CM = self.phaseSpace == "CM"
		return retVal

	def checkPole(self, position = None, delta = 1.e-5, verbose = False):
		fcn             = self.produceFunction()
		fcn.secondSheet = True
		if position is None:
			try:
				polePar = self.mRho**2 + 1.j*self.mRho*self.Grho
			except Exception as e:
				print "noPar:",self.inFileName
				return False
		else:
			polePar = position
		rePlus          = fcn.complexCall(polePar + delta)
		reMinus         = fcn.complexCall(polePar - delta)
		retVal          = True
		if rePlus.real*reMinus.real > 0.:
			if verbose:
				print "No sign change at pole position: Real part in real direction"
			retVal  = False
		if rePlus.imag*reMinus.imag > 0.:
			if verbose:
				print "No sign change at pole position: Imag part in real direction"
			retVal  = False
		imPlus          = fcn.complexCall(polePar + 1.j*delta)
		imMinus         = fcn.complexCall(polePar - 1.j*delta)
		if imPlus.real*imMinus.real > 0.:
			if verbose:
				print "No sign change at pole position: Real part in imag direction"
			retVal  = False
		if imPlus.imag*imMinus.imag > 0.:
			if verbose:
				print "No sign change at pole position: Imag part in imag direction"
			retVal  = False
		if retVal:
			if verbose:
				print "Pole is good: ", rePlus,reMinus,imPlus,imMinus
		else:
			if verbose:
				print "Pole is bad:  ", rePlus,reMinus,imPlus,imMinus
		return retVal

	def circleIntegral(self, position, radius, nPoints):
		step = 2*pi/nPoints
		val  = 0.
		fcn  = self.produceFunction()
		fcn.secondSheet = True
		for i in range(nPoints):
			z    = position + radius*(cos(i*step) + 1.j*sin(i*step))
			dz   = step*radius*(-sin(i+step) + 1.j *cos(i*step))
			val += fcn.complexCall(z)*dz
		return val

	def boxBelowRealAxis(self, zMin, zMax, nPoints, epsilon = 1.e-6):
		halfLength      = (zMax-zMin)/2
		dz              = (zMax-zMin)/nPoints
		val             = 0.
		fcn             = self.produceFunction()
		fcn.secondSheet = True
		for i in range(nPoints):
			re   = zMin + (i+.5)*dz
			z0   = re - 1.j*epsilon
			z1   = re - 1.j*halfLength
			val += (fcn.complexCall(z1) - fcn.complexCall(z0))*dz
		for i in range(nPoints/2):
			im   = (i+.5)*dz
			z0   = zMin + im
			z1   = zMax + im
			val += (fcn.complexCall(z1) - fcn.complexCall(z0))*1.j*dz
		return val

class resultHolder:
	def __init__(self, folder, numLim = 1.e-13, fileNamePreFilter = []):
		self.folder  = folder
		self.filters = []
		self.numLim  = numLim
		self.load(fileNamePreFilter = fileNamePreFilter)

	def load(self, verbose = False, fileNamePreFilter = []):
		self.results = []
		for fn in os.listdir(self.folder):
			reject = False
			for fnpf in fileNamePreFilter:
				if not fnpf in fn:
					reject = True
					break
			if reject:
				continue
			if "Kmatrix" in fn and fn.endswith(".dat"):
				try:
					self.results.append(KmatrixFitResult(self.folder+os.sep+fn, self.numLim))
					self.results[-1].produceFunction()
				except ValueError as e:
					if verbose:
						print str(e)
					continue

	def getFilteredList(self):
		retVal = []
		for result in self.results:
			good = True
			for f in self.filters:
				if not f(result):
					good = False
					break
			if good:
				retVal.append(result)
		return retVal

def filterCM(result):
	return result.phaseSpace == "CM"

def filterRho(result):
	return result.phaseSpace == "rho"

def filter1pp(result):
	return result.wave == "1pp"

def filter1mp(result):
	return result.wave == "1mp"

def filterRightConfiguration(result):
	if not result.kPol == 7:
		return False
	if not result.pPOl == (7,4):
		return False
	return True

class filterRange2:
	def __init__(self,val):
		self.val = val

	def __call__(self, result):
		return result.range == self.val

class filterT:
	def __init__(self,tBin):
		self.t = tBin

	def __call__(self, result):
		return result.tBin == self.t
				
class filterRange3:
	def __init__(self, mMin = None, mMax = None):
		self.min = mMin
		self.max = mMax
	
	def __call__(self, result):
		if self.min is not None:
			if not self.min == result.mRange[0]:
				return False
		if self.max is not None:
			return self.max == result.mRange[1]	
		print "WARNING: filterMrange(None, None) has no effect"
		return True

def doInfo(vals, canv, minDiff = 1.e-4):
	if len(vals) < 2:
		return False
	mMin =  float('inf')			
	mMax = -float('inf')
	Gmin =  float('inf')
	Gmax = -float('inf')
	for v in vals:
		mMin = min(mMin, v[1])
		mMax = max(mMax, v[1])
		Gmin = min(Gmin, v[2])
		Gmax = max(Gmax, v[2])
	mRange = mMax-mMin
	Grange = Gmax-Gmin
	nBins  = 100
	margin = .1
	hist   = ROOT.TH2D("hol","hol", nBins, mMin-margin*mRange, mMax+margin*mRange, nBins, Gmin-margin*Grange, Gmax+margin*Grange)
	bestMbin = hist.GetXaxis().FindBin(vals[0][1])
	mBinUp   = hist.GetXaxis().FindBin(vals[0][1]+.005)
	mBinLow  = hist.GetXaxis().FindBin(vals[0][1]-.005)

	bestGbin = hist.GetYaxis().FindBin(vals[0][2])
	GbinUp   = hist.GetYaxis().FindBin(vals[0][2]+.005)
	GbinLow  = hist.GetYaxis().FindBin(vals[0][2]-.005)

	for v in vals:
		print mMin, mMax, v[1]
		print Gmin, Gmax, v[2]
		hist.Fill(v[1], v[2])
	for i in range(nBins):
		if hist.GetBinContent(bestMbin, i+1) == 0.:
			hist.SetBinContent(bestMbin, i+1, .1)
			hist.SetBinContent(mBinUp, i+1, .05)
			hist.SetBinContent(mBinLow, i+1, .05)
		if hist.GetBinContent(i+1, bestGbin) == 0.:
			hist.SetBinContent(i+1, bestGbin, .1)
			hist.SetBinContent(i+1, GbinUp, .05)
			hist.SetBinContent(i+1, GbinLow, .05)
	hist.Draw("COL")
	canv.Update()
	raw_input()

def no13(result):
	return not result.mRange[0] == 13

class rhoPrimeRangeFilter:
	def __init__(self, mRange, Grange):
		self.mRange = mRange
		self.Grange = Grange

	def __call__(self, res):
		if hasattr(res, "mRhoPrime"):
			if res.mRhoPrime < self.mRange[0]:
				return False
			if res.mRhoPrime > self.mRange[1]:
				return False
		else:
			return False
		if hasattr(res, "GrhoPrime"):
			if res.GrhoPrime < self.Grange[0]:
				return False
			if res.GrhoPrime > self.Grange[1]:
				return False
		else:
			return False
		return True

def main():

# # #  Filters

	ROOT.gStyle.SetOptStat(0)
	twoPoleFilter      = lambda x: x.nPol == 2
	mRange1Filter      = lambda x: x.mRange[0] + 1 == x.mRange[1]
	filterComa         = lambda x: x.coma     is not None
	hasBW              = lambda x: x.BWresult is not None
	filterM3range      = filterRange3(14,50)
	range2Filter       = filterRange2(None)
	tFilter            = filterT(0)
	rprf               = rhoPrimeRangeFilter((1.1,3.),(.1,.7))
# # # # # All results

	allResults         = resultHolder("./KmatrixResults", 1.e-11, fileNamePreFilter = ["1mp_"])
	allResults.filters = [filter1mp, hasBW, filterCM]
	fl                 = allResults.getFilteredList()

	mHist = ROOT.TH1D("rhoMass", "rhoMass",  300, 0.5, 1. )
	Ghist = ROOT.TH1D("rhoWidth","rhoWidth", 300, 0. ,  .3)

	mBWhist = ROOT.TH1D("rhoBWmass", "rhoBWmass",  300, 0.5, 1. )
	GbwHist = ROOT.TH1D("rhoBWwidth","rhoBWwidth", 300, 0. ,  .3)

	mBWhist.SetLineColor(2)
	GbwHist.SetLineColor(2)
	matcessible = {}
	for t in range(4):
		for m in range(50):
			matcessible[(t,m)] = []
	for f in fl:
		matcessible[(f.tBin, f.mRange[0])].append((f.chi2, f.mRho, f.Grho, f.BWresult[0][0],f.BWresult[1][0]))
	for bin in matcessible:
		matcessible[bin].sort()
#	c = ROOT.TCanvas()
#	for t in range(4):
#		for m in range(50):
#			doInfo(matcessible[(t,m)],c)
#	return

	for bin in matcessible:
		if len(matcessible[bin]) == 0:
			continue
		mHist.Fill(matcessible[bin][0][1])
		Ghist.Fill(matcessible[bin][0][2])
		mBWhist.Fill(matcessible[bin][0][3])
		GbwHist.Fill(matcessible[bin][0][4])

#	for f in fl:
#		for g in fl:
#			if f.mRange[0] == g.mRange[0] and f.tBin == g.tBin:
#				print "========================================"
#				print f.mRho-g.mRho, f.Grho-g.Grho
#				print f.BWresult[0][0]-g.BWresult[0][0],f.BWresult[1][0]-g.BWresult[1][0]
#		mGcoma = f.makeMassWidthComa()
#		mHist.Fill(f.mRho)
#		Ghist.Fill(f.Grho)
#
#		mBWhist.Fill(f.BWresult[0][0])
#		GbwHist.Fill(f.BWresult[1][0])

#		print f.mRho,mGcoma[0][0]**.5,f.Grho,mGcoma[1][1]**.5, f.inFileName
#		print f.mRho**2, f.coma[0][0], -f.Grho*f.mRho, f.coma[1][1]

	mHist.Draw()
	mBWhist.Draw("SAME")
	raw_input()
	Ghist.Draw()
	GbwHist.Draw("SAME")
	raw_input()

	return

	nPrime             = 0
	c2s                = []
	histRhoPrime       = ROOT.TH2D("rp", "rp", 100, 0., 2., 100, 0., 1.6)
	for res in fl:

#		print "-------------------------------------------------------------------"
#		print res.inFileName
#		print res.chi2/res.ndf,"chi2"
		if hasattr(res, 'fValPrime'):
			if abs(res.fValPrime) < 1.e-11:
				c2s.append(res.chi2/res.ndf)
				print "-------------------------------------------------------------------"

				primePos = res.mRhoPrime**2 + 1.j*res.mRhoPrime*res.GrhoPrime
				print "(m/G)'",res.mRhoPrime,res.GrhoPrime,res.fValPrime/res.ndf
				histRhoPrime.Fill(res.mRhoPrime,res.GrhoPrime)
				print res.inFileName, res.chi2
		
#				print res.mRhoPrime, res.GrhoPrime, "< < < <", res.fValPrime
				if not res.checkPole(primePos, verbose = False, delta = 1.e-3):
					pass
	#				func = res.produceFunction()
	#				func.secondSheet = True
	#				dd   = 0.1
	#				nn   = 200
	#				step = 2*dd/nn
	#				histReRe = ROOT.TH1D('rere','rere',nn,primePos.real-dd, primePos.real+dd)
	#				histReIm = ROOT.TH1D('rere','rere',nn,primePos.real-dd, primePos.real+dd)
	#				histImRe = ROOT.TH1D('rere','rere',nn,primePos.real-dd, primePos.real+dd)
	#				histImIm = ROOT.TH1D('rere','rere',nn,primePos.real-dd, primePos.real+dd)
	#				for i in range(nn):
	#					valRe = func.complexCall(primePos -      dd + (i+0.5)*step)
	#					valIm = func.complexCall(primePos - 1.j*(dd + (i+0.5)*step))
	#					histReRe.SetBinContent(i+1, valRe.real)
	#					histReIm.SetBinContent(i+1, valRe.imag)
	#					histImRe.SetBinContent(i+1, valRe.real)
	#					histImIm.SetBinContent(i+1, valRe.imag)
	#				histReRe.Draw()
	#				raw_input()
	#				histReIm.Draw()
	#				raw_input()
	#				histImRe.Draw()
	#				raw_input()
	#				histImIm.Draw()
	#				raw_input()
				else:
					primePos = res.mRhoPrime**2 + 1.j*res.mRhoPrime*res.GrhoPrime
					nPrime += 1

#		print
	print len(fl),nPrime
	c2s.sort()
	chi2hist = ROOT.TH1D("c2","c2", len(c2s), 0., c2s[-1])
	for c2 in c2s:
		chi2hist.Fill(c2)
	histRhoPrime.Draw("COLZ")
	raw_input()
	chi2hist.Draw()
	raw_input()

if __name__ == "__main__":
	main()
