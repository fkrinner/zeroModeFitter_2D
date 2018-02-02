import numpy as np
import numpy.linalg as la
import os, sys 
import pyRootPwa

from math import pi, exp, log, atan

numLim = 1.E-10
INF = float("inf")

def getNDp(sector):
	"""
	Gets the number of freed isobar bins for D decays
	"""
	nDpfile = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/build/nDp.dat"
	ns = []
	with open(nDpfile, 'r') as inin:
		for line in inin.readlines():
			ns += [int(chunk) for chunk in line.split()]
	if not len(ns) == 3:
		raise IOError("Dit not get exactly 3 values")
	if "0++" in sector:
		print "HARDCODE::::: getNdp!!!!!!!!!!!!!!"
		return 49
		return ns[0]
	elif "1--" in sector:
		print "HARDCODE::::: getNdp!!!!!!!!!!!!!!"
		return 43
		return ns[1]
	elif "2++" in sector:
		return ns[2]
	else:
		raise ValueError("Unknown sector '" + sector + "'")	

nF0   = 62
nRho  = 56
nF2   = 56
nRho3 = 56

def zeroForSectors(sectors, title):
	"""
	Determines, whether a zero mode contributes to the sector at hand
	"""
	for sector in sectors:
		if sector in title:
			return True
	return False

def countSectors(hist):
	"""
	Counts the sectors, a zero mode is contributung to
	"""
	return len(getZeroHistSectors(hist))

def getZeroHistSectors(hist):
	"""
	Gets the list of sectors for a zero mode
	"""
	return  hist.GetTitle().split("<|sec++>")[0].split("<|>")

def getNforSector(sector):
	"""
	Gets the number of bind for a sctor
	"""
	if "Dp" in sector:
		return getNDp(sector)
	elif "[pi,pi]0++" in sector:
		return nF0
	elif "[pi,pi]1--" in sector:
		return nRho
	elif "[pi,pi]2++" in sector:
		return nF2
	elif "[pi,pi]3--" in sector:
		return nRho3
	else:
		raise ValueError("Could nto detect nBins for '" + sector + "'")

def getZeroHistBorders(zeroHist):
	"""
	Gets the sector borders for a zero mode histogram
	"""
	sectors = getZeroHistSectors(zeroHist)
	borders = [0]
	for sector in sectors:
		borders.append(borders[-1] + getNforSector(sector))
	return borders

def renormToBinWidth(hist, exponent = 1.):
	"""
	Renorms a histogram to the bin width (with a certain exponent)
	"""
	for i in range(hist.GetNbinsY()):
		binWidth = hist.GetYaxis().GetBinWidth(i+1)
		for j in range(hist.GetNbinsX()):
			hist.SetBinContent(j+1,i+1, hist.GetBinContent(j+1,i+1)/binWidth**exponent)
			hist.SetBinError(j+1, i+1, hist.GetBinError(j+1, i+1)/binWidth**exponent)

def pinv(matrix, numLim = 1.e-13):
	"""
	Own method fot the pseudo inverse of a matrix
	"""	
	dim = len(matrix)
	val, vec = la.eig(matrix)
	for i in range(dim):
		if abs(val[i]) < numLim:
			val[i] = 0.
		elif val[i].real < -2.e-8:
			raise ValueError("Negative eingenvalue: " + str(val[i]))
		else:
			val[i] = 1./val[i]
	return np.dot(vec, np.dot(np.diag(val), np.transpose(vec)))

def isValidPhaseSpace(m3pi, m2pi, mPi = 0.139):
	"""
	Ckecks, if the process is kinametically valid
	"""
	if m3pi >= m2pi + mPi:
		return True
	return False

def LtoInt(L):
	"""
	Gets the spin quantum number for the corresponding letter
	"""
	if L == 'S':
		return 0
	if L == 'P':
		return 1
	if L == 'D':
		return 2
	if L == 'F':
		return 3
	if L == 'G':
		return 4
	if L == 'H':
		return 5


def getOnePhaseDirection(ampl):
	"""
	Returns the direction of a phase roatation of the complex inpit number
	"""
	return (-ampl.imag, ampl.real)

def normVector(vector):
	"""
	Normalizes a the vector to 1
	"""
	dim  = len(vector)
	norm = 0.
	for v in vector:
		norm += v**2
	norm**=.5
	for i in range(dim):
		vector[i] /= norm
	return vector

def getPhaseDirection(ampls):
	"""
	Gets the normalized direction of a global phase rotation in all complex numbers
	"""
	dim = len(ampls)
	if dim%2 != 0:
		raise ValueError("Odd number of values... Error, cannot determine phase-direction")
	retVal = np.zeros(dim)
	for i in range(dim/2):
		pd     = getOnePhaseDirection(ampls[2*i] + 1.j * ampls[2*i+1])
		retVal[2*i  ] = pd[0]
		retVal[2*i+1] = pd[1]
	return normVector(retVal)

def getZeroModeNumber(hist):
	"""
	Returns the ID number of a zero mode histograms
	"""
	name = hist.GetName()
	if not name.startswith('zero'):
		raise NameError("'" + name + "' is not the name of a zeroMode histogram")
	return int(name.split('_')[0][4:])

def sumUp(ll):
	"""
	Sums up the scalar content of an n-dimensional array
	"""
	if hasattr(ll, '__len__'):
		summ = 0.
		for l in ll:
			summ += sumUp(l)
		return summ
	return ll

def cloneZeros(lst):
	"""
	Clones the shape of an array, filled with zeros
	"""
	retVal = []
	for val in lst:
		if hasattr(val, '__len__'):
			retVal.append(cloneZeros(val))
		else:
			retVal.append(0.)
	return retVal

def addBtoA(A,B, weight = 1., noBelowZero = False):
	"""
	Adds the entries od array B at the corresponding palces of the equal shaped array A, unsing the weight as factor
	"""
	if not len(A) == len(B):
		raise ValueError("addBtoA(...): Size mismatch ('" + str(len(A)) + "' != '" + str(len(B)) + "')" )
	for i in range(len(A)):
		if hasattr(A[i], '__len__'):
			addBtoA(A[i],B[i], weight, noBelowZero = noBelowZero)
		else:
			if noBelowZero and B[i] < 0.:
				continue
			A[i] += weight*B[i]

def divideAbyB(A,B, ignoreZero = False):
	"""
	Divides the entries of A by the corresponding entris of B
	"""
	if not len(A) == len(B):
		raise ValueError("divideAbyB(...): Size mismatch ('" + str(len(A)) + "' != '" + str(len(B)) + "')" )
	for i in range(len(A)):
		if hasattr(A[i], "__len__"):
			divideAbyB(A[i], B[i], ignoreZero = ignoreZero)
		else:
			if ignoreZero and B[i] == 0.:
				A[i] == 0.
			else:
				A[i] /= B[i]

def weightedSum(weights, params):
	"""
	Calculates the weightd sum of the arrays
	"""
	if len(params) == 0:
		return []
	for m in weights:
		retVal = cloneZeros(params[m])
		break
	for m in weights:
		addBtoA(retVal, params[m], weights[m])
	return retVal

def invertAllEntries(A, nanToZero = False):
	for i in range(len(A)):
		if hasattr(A[i], "__len__"):
			invertAllEntries(A[i], nanToZero = nanToZero)
		else:
			if nanToZero and (A[i] == 0. or np.isnan(A[i])):
				A[i] = 0.
			else:
				A[i] = 1./A[i]

def printZeroStructure(matrix):
	dim = len(matrix)
	for i in range(dim):
		st = ""
		for j in range(dim):
			if matrix[i,j] == 0.:
				st += '0'
			else:
				st += '1'
		print st

def checkLaTeX():
	"""
	Checks, if 'uplatex' has been called. Raises an exception, if not, 
	so the program terminates right there and does not keep calculating 
	and then crash on plotting.
	"""
	if not "/nfs/mnemosyne/sys/slc6/contrib/texlive/2013/bin/x86_64-linux" in os.environ["PATH"]:
		raise RuntimeError("LaTeX not set correctly, aborting")

globalHistCount = 0
def get3PiHistogram(name = ""):
	if name == "":
		global globalHistCount
		name = "hist_"+str(globalHistCount)
		globalHistCount += 1
	hist = pyRootPwa.ROOT.TH1D(name, name, 50, 0.5,2.5)
	hist.SetDirectory(0)
	return hist

def loadAmplsTM(inFileName):
	"""
	Loads amplitudes from a text file with lines: tBin mBin real imag
	Assumes four t' bins and 50 m bins
	"""
	ampls = [[0. for _ in range(50)] for __ in range(4)]
	with open(inFileName, 'r') as inin:
		count = 0
		for line in inin.readlines():
			chunks = line.split()
			tBin   = int(chunks[0])
			mBin   = int(chunks[1])
			real   = float(chunks[2])
			imag   = float(chunks[3])
			if not ampls[tBin][mBin] == 0.:
				raise ValueError("Amplitude for "+str(tBin)+" "+str(mBin)+"set twice")
			ampls[tBin][mBin] = real + 1.j*imag
			count += 1
	if not count == 200:
		raise ValueError("Not all amplitudes set")
	return ampls

def changeReferenceWave(histListReal, histListImag, histListIndx, comaHists, refHistReal, refHistImag, refHistIndex, startBin, stopBin):
	"""
	Transforms the COMA to another reference wave (But only for the indices in the given histograms)
	"""
	nHists = len(histListReal)
	for iComa, iHist in enumerate(range(startBin, stopBin)):
		dim      = comaHists[iComa].GetNbinsX()	
		vals     = np.zeros(dim)
		jacobian = np.zeros((dim,dim))
		coma     = np.zeros((dim,dim))
		for i in range(dim):
			for j in range(dim):
				coma[i,j] = comaHists[iComa].GetBinContent(i+1, j+1)
		refIndex = refHistIndex.GetBinContent(iHist+1)
		reRef  = refHistReal.GetBinContent(iHist+1)
		imRef  = refHistImag.GetBinContent(iHist+1)
		absRef = (reRef**2+imRef**2)**.5
		vals[2*refIndex  ] = reRef
		vals[2*refIndex+1] = imRef
		for h in range(len(histListReal)):
			for i in range(histListReal[h].GetNbinsY()):
				index = histListIndx[h].GetBinContent(iHist+1, i+1)
				if index == 0:
					break
				vals[2*index  ] = histListReal[h].GetBinContent(iHist+1, i+1)
				vals[2*index+1] = histListImag[h].GetBinContent(iHist+1, i+1)
				compl = (histListReal[h].GetBinContent(iHist+1, i+1) + 1.j *histListImag[h].GetBinContent(iHist+1, i+1))*(reRef-1.j*imRef)/absRef
				histListReal[h].SetBinContent(iHist+1, i+1, compl.real)
				histListImag[h].SetBinContent(iHist+1, i+1, compl.imag)



		for i in range(dim/2): # *= (reRef - i imRef)/absRef
			jacobian[2*i  ,2*i  ] = reRef/absRef # dReNew/dReOld
			jacobian[2*i  ,2*i+1] = imRef/absRef # dReNew/dImOld
			jacobian[2*i+1,2*i  ] =-imRef/absRef # dImNew/dReOld
			jacobian[2*i+1,2*i+1] = reRef/absRef # dImNew/dImOld
			
			jacobian[2*i  ,2*refIndex  ] = vals[2*i  ]/absRef - (vals[2*i  ]*reRef + vals[2*i+1]*imRef)*reRef/absRef**3
			jacobian[2*i  ,2*refIndex+1] = vals[2*i+1]/absRef - (vals[2*i  ]*reRef + vals[2*i+1]*imRef)*imRef/absRef**3
			jacobian[2*i+1,2*refIndex  ] = vals[2*i+1]/absRef - (vals[2*i+1]*reRef - vals[2*i  ]*imRef)*reRef/absRef**3
			jacobian[2*i+1,2*refIndex+1] =-vals[2*i  ]/absRef - (vals[2*i+1]*reRef - vals[2*i  ]*imRef)*imRef/absRef**3
	

		intermed = np.dot(jacobian, coma)
		newComa  = np.dot(intermed, np.transpose(jacobian))

		for i in range(dim):
			for j in range(dim):
				comaHists[iComa].SetBinContent(i+1,j+1, newComa[i,j])

def gaus1D(x,par):
	"""
	Parameter ordering: norm, meanX, sigX
	"""
	exponent = -(x-par[1])**2/par[2]**2/2
	norm     =  par[0]/(2*pi)**.5/par[2]
	retVal   =  norm*exp(exponent)
	return retVal

def gaus2D(x,y,par):
	"""
	Parameter ordering: norm, meanX, sigX, meanY, sigY, correl
	"""
	if not len(par) == 6:
		raise ValueError("Wrong number of parameters: "+str(len(par)))
	correl = atan(par[5])*2/pi*par[2]*par[4]
#	correl   = par[5]*par[2]*par[4]
	detComa  = par[2]**2 * par[4]**2 - correl**2
	inverse  = [[par[4]**2/detComa, -correl/detComa],[-correl/detComa,par[2]**2/detComa]]
	exponent = 0.
	d        = [par[1] - x, par[3] - y]
	for i in range(2):
		for j in range(2):
			exponent += d[i]*d[j]*inverse[i][j]
	retVal = par[0]/2/pi/detComa**.5 * exp(-exponent/2)
#	if isnan(retVal):
#		print "NAN",  detComa, exponent, retVal
	return retVal

def gaus3D(x,y,z,par):
	"""
	Parameter ordering: norm, meanX, sigX, meanY, sigY, meanZ, sigZ, correlXY, correlXZ, correlYZ
	"""
	if not len(par) == 10:
		raise ValueError("Wrong number of parameters: "+str(len(par)))
	cXY = atan(par[7])*2/pi*par[2]*par[4]
	cXZ = atan(par[8])*2/pi*par[2]*par[6]
	cYZ = atan(par[9])*2/pi*par[4]*par[6]
	d        = [par[1] - x, par[3] - y, par[5] - z]
	coma = np.asarray([[par[2]**2, cXY, cXZ],[cXY, par[4]**2, cYZ],[cXZ, cYZ, par[6]**2]])
	inverse = la.inv(coma)
	exponent = 0.
	for i in range(3):
		for j in range(3):
			exponent += d[i]*d[j]*inverse[i,j]
	retVal = par[0]/(8*(pi)**3*la.det(coma))**.5 * exp(-exponent/2)
	return retVal
	

class oneDgaussFuncLike:
	def __init__(self, values):
		self.values  = values
		self.nPoints = len(values)
		if self.nPoints == 0:
			raise ValueError("No values given")
		if isinstance(values[0], float):
			self.isFloat  = True
			self.weighted = False
		else:
			self.isFloat = False
			if len(values[0]) == 1:
				self.weighted = False
			elif len(values[0]) == 2:
				self.weighted = True
			else:
				raise ValueError("Number of values invalid: Neither x, (x), nor (x,weight)")
		self.nPrint = -1
		self.count  =  0

	def __call__(self,par):
		like = 0.
		self.count += 1
		pars = [1.]*3
		for i in range(2):
			pars[i+1] = par[i]
		for vals in self.values:
			if self.isFloat:
				x = vals
			else:
				x = vals[0]
			gausVal = gaus1D(x, pars)
			if gausVal < 0.:
				continue
			logVal = log(gausVal)
			if self.weighted:
				like -= vals[1] * logVal
			else:
				like -= logVal
		if self.nPrint > 0 and self.count%self.nPrint == 0:
			print '#' + str(self.count), 
			print like,
			print par
		return like			

class twoDgaussFuncLike:
	def __init__(self, values):
		self.values  = values
		self.nPoints = len(values)
		if self.nPoints == 0:
			raise ValueError("No values given")
		if len(self.values[0]) == 2:
			self.weighted = False
		elif len(self.values[0]) == 3:
			self.weighted = True
		else:
			raise ValueError("Number of values invalid: Neither (x,y), nor (x,y,weight)")
		self.nPrint = -1
		self.count  =  0

	def __call__(self, par):
		like = 0.
		self.count += 1
		pars = [1.]*6
		for i in range(5):
			pars[i+1] = par[i]
		for vals in self.values:
			gausVal = gaus2D(vals[0], vals[1], pars)
			if gausVal <= 0.:
				continue
			logVal = log(gausVal)
			if self.weighted:
				like -= vals[2] * logVal
			else:
				like -= logVal
		if self.nPrint > 0 and self.count%self.nPrint == 0:
			print '#' + str(self.count), 
			print like,
			print par
		return like

class threeDgaussFuncLike:
	def __init__(self, values):
		self.values  = values
		self.nPoints = len(values)
		if self.nPoints == 0:
			raise ValueError("No values given")
		if len(self.values[0]) == 3:
			self.weighted = False
		elif len(self.values[0]) == 4:
			self.weighted = True
		else:
			raise ValueError("Number of values invalid: Neither (x,y), nor (x,y,weight)")
		self.nPrint = -1
		self.count  =  0

	def __call__(self, par):
		like = 0.
		self.count += 1
		pars = [1.]*10
		for i in range(9):
			pars[i+1] = par[i]
		for vals in self.values:
			gausVal = gaus3D(vals[0], vals[1],vals[2], pars)
			if gausVal <= 0.:
				continue
			logVal = log(gausVal)
			if self.weighted:
				like -= vals[3] * logVal
			else:
				like -= logVal
		if self.nPrint > 0 and self.count%self.nPrint == 0:
			print '#' + str(self.count), 
			print like,
			print par
		return like

class twoDgausFunc:
	def __init__(self, hist2D, isLike = False):
		self.hist   = hist2D
		self.isLike = isLike
		self.count  =  0
		self.nPrint = -1
		self.mbc    = -1 # minBinContent for evaluation

	def __call__(self, par):
		chi2 = 0.
		like = 0.
		pars = par
		self.count += 1
		if self.isLike:
			pars = [1.]*6
			for i in range(5):
				pars[i+1] = par[i]

		for i in range(self.hist.GetNbinsX()):
			x = self.hist.GetXaxis().GetBinCenter(i+1)
			for j in range(self.hist.GetNbinsY()):
				entries = self.hist.GetBinContent(i+1, j+1)
				if entries == 0. and self.isLike:
					continue
				if entries < self.mbc:
					continue

				y = self.hist.GetYaxis().GetBinCenter(j+1)
				gausVal = gaus2D(x,y,pars)
				if not entries == 0.:
					chi2 += (entries - gausVal)**2/entries
				else:
					chi2 += gausVal**2
					pass
				if gausVal > 0.:
					like -= entries*log(gausVal)
		if self.nPrint > 0 and self.count%self.nPrint == 0:
			print '#' + str(self.count), 
			if self.isLike:
				print like,
			else:
				print chi2
			print par

		if self.isLike:
			return like
		else:
			return chi2

	def ndf(self):
		ndf = 0
		for i in range(self.hist.GetNbinsX()):
			for j in range(self.hist.GetNbinsY()):
				if self.hist.GetBinContent(i+1,j+1) >= self.mbc:
					ndf += 1
		return ndf - 6

	def makeDrawHist(self, par):
		pars = par
		if self.isLike:
			pars = [1.]*6
			for i in range(5):
				pars[i+1] = par[i]
		newHist = self.hist.Clone()
		for i in range(self.hist.GetNbinsX()):
			x = self.hist.GetXaxis().GetBinCenter(i+1)
			for j in range(self.hist.GetNbinsY()):
				y = self.hist.GetYaxis().GetBinCenter(j+1)
				gausVal = gaus2D(x,y,pars)
				newHist.SetBinContent(i+1, j+1, gausVal)
		return newHist


def main():
	from random import random
	def leDiste():
		while True:
			x = (random() - .5)
			if exp(-abs(x)) > random():
				return x

	nPoints = 10000
	pts     = []
	avg = 0.
	for i in range(nPoints):
		print i
		p = leDiste()
		pts.append([p,random()])
		avg += p
	avg /= len(pts)
	sig = 0.
	for p in pts:
		sig += (avg - p[0])**2
	sig/= len(pts)-1
	sig**=.5
	print "avg",avg,"sig",sig

	likeFunc = oneDgaussFuncLike(pts)
	par = [0.,1.]
	import scipy.optimize as opt
	res = opt.minimize(likeFunc, par)
	print res.x


if __name__ == "__main__":
	main()

	

