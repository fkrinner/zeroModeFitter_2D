#!/usr/bin/python
# vectorFormFactor.py
# Created: 2018-07-16 11:16:22.825312
# Author: Fabian Krinner
import os, sys
from cmath import exp, pi
from cmath import log as clog
from math import atan, atan2, log
#from globalDefinitions import mPi, mK, mTau
from globalDefinitions import mTau

from convertedFortranCode import Ftilde

from parameterTrackingParameterizations import parametersTrackingParameterization
import scipy.integrate as integrate
import numpy as np

 # # # # # THIS PARAMETER ORDERING IS FIXED NOW... NEVER EVER CHANGE # # # !!!
#	lambPrime  = param[ 0]
#	lambPPrime = param[ 1]
 # # # # # THIS PARAMETER ORDERING IS FIXED NOW... NEVER EVER CHANGE # # # !!!
#	Mrho       = param[ 2]
#	Grho       = param[ 3]
#	MrhoPrime  = param[ 4]
#	GrhoPrime  = param[ 5]
#	MrhoPPrime = param[ 6]
#	GrhoPPrime = param[ 7]
 # # # # # THIS PARAMETER ORDERING IS FIXED NOW... NEVER EVER CHANGE # # # !!!
#	alphPrime  = param[ 8]
#	phiPrime   = param[ 9]
#	aplhPPrime = param[10]
#	phiPPrime  = param[11]
 # # # # # THIS PARAMETER ORDERING IS FIXED NOW... NEVER EVER CHANGE # # # !!!

mK  = 0.4957
mPi = 0.13957018

def conformal(x,x0):
	return x0*atan(x/x0)

def Fv(s, param, renScale2 = .77**2, sCut = mTau**2, cutVal = pi, x0 = 10., allowSubThr = False):
	mPi2     = mPi**2
	exponent =  param[0]*s/mPi2 + (param[1] - param[0]**2)*(s/mPi2)**2/2 + s**3/pi * dispersiveIntegral(s, param, renScale2 = renScale2, sCut = sCut, cutVal = cutVal, allowSubThr = allowSubThr)
#	print s,'|',param[0]*s/mPi2,  (param[1] - param[1]**2)*(s/mPi2)**2/2,  s**3/pi * dispersiveIntegral(s, param, renScale2 = renScale2, sCut = sCut, cutVal = cutVal, allowSubThr = allowSubThr).real,'|',exponent.real
	if x0 is None:
		return exp(exponent)
	else:
#@		print exponent, conformal(exponent.real, x0) + 1.j*exponent.imag, exp(conformal(exponent.real, x0) + 1.j*exponent.imag), " }}}}}}}}}}}}}}{{}}}}}}}}}}}}}}}}}}}}}{{"
		return exp(conformal(exponent.real, x0) + 1.j*exponent.imag)

def FvTilde(s, param, renScale2, allowSubThr = False):

	part1 = (param[2]**2 + (param[8]*exp(1.j*param[9]) + param[10]*exp(1.j*param[11]))*s)/(param[2]**2 - s+kappaRho(s,param, allowSubThr = allowSubThr) * (A(s, mPi, renScale2, allowSubThr = allowSubThr) + A(s, mK, renScale2, allowSubThr = allowSubThr)/2))
	part2 = -param[ 8]*exp(1.j*param[ 9])*s/Den(s, param, 1, renScale2, allowSubThr = allowSubThr)
	part3 = -param[10]*exp(1.j*param[11])*s/Den(s, param, 2, renScale2, allowSubThr = allowSubThr)
	return part1 + part2 + part3

def phiV(s,param, renScale2, allowSubThr = False):
	ampl = FvTilde(s, param, renScale2, allowSubThr = allowSubThr)
#	return atan2(ampl.imag, ampl.real)
	phi = atan(ampl.imag/ampl.real)
	if ampl.real < 0. or ampl.imag < 0.:
		phi += pi
	return phi

def sigma(s, mPart, allowSubThr = False):
	xpr = 1. - 4*mPart**2/s + 0.j
	if xpr.real < 0. and not allowSubThr:
		return 0.
	return xpr**.5

def A(s, mPart, renScale2, allowSubThr = False):
	sig    = sigma(s, mPart, allowSubThr = allowSubThr)
	retVal = log(mPart**2/renScale2) + 8*mPart**2/s - 5./3. + sig**3*clog((sig+1)/(sig-1.)) # +1.j pi to obtain the real part (at least above threshold)
	return retVal

#def Areal(s, mPart, renScale2, allowSubThr = False):
#	sig    = sigma(s, mPart, allowSubThr = allowSubThr)
#	print s, mPart,sig
#	retVal = log(mPart**2/renScale2) + 8*mPart**2/s - 5./3. + sig**3*(clog((sig+1)/(sig-1.)) + 1.j*pi) # +1.j pi to obtain the real part (at least above threshold)
#	print retVal.imag, pi*sig**3
#	return retVal

def Den(s, param, nPrime, renScale2, allowSubThr = False):
#	retVal = param[2+2*nPrime]**2 - s + kappaPrime(s, param, nPrime, allowSubThr = allowSubThr) *Areal(s, mPi, renScale2, allowSubThr = allowSubThr) - 1.j*param[2+2*nPrime]*gammaPrime(s, param, nPrime, allowSubThr = allowSubThr)
	retVal = param[2+2*nPrime]**2 - s + kappaPrime(s, param, nPrime, allowSubThr = allowSubThr) *A(s, mPi, renScale2, allowSubThr = allowSubThr)
	return retVal

def kappaRho(s,param, allowSubThr = False):
	sig = sigma(param[2]**2, mPi, allowSubThr = allowSubThr)**3 + sigma(param[2]**2, mK, allowSubThr = allowSubThr)**3/2
	if sig == 0.:
		return 0.
	retVal = param[3]/param[2]*s/(pi*sig)
	return retVal

def kappaPrime(s,param,nPrime, allowSubThr = False):
	sig = sigma(param[2+2*nPrime]**2, mPi, allowSubThr = allowSubThr)
	if sig == 0.:
		return 0.
	retVal = param[3+2*nPrime]/param[2+2*nPrime] *s/(pi*sig**3)
	return retVal

#def gammaRho(s,param, allowSubThr = False):
#	sig = sigma(param[2]**2,mPi, allowSubThr = allowSubThr)**3 + sigma(param[2]**2,mK, allowSubThr = allowSubThr)**3/2
#	if sig == 0.:
#		return 0.
#	retVal = param[3]*s/param[2]**2 *(sigma(s,mPi, allowSubThr = allowSubThr)**3 + sigma(s,mK, allowSubThr = allowSubThr)**3/2)/sig
#	return retVal

#def gammaPrime(s, param, nPrime, allowSubThr = False):
#	sig = sigma(param[2+2*nPrime]**2, mPi, allowSubThr = allowSubThr)
#	if sig == 0.:
#		return 0.
#	retVal = param[3+2*nPrime]*s/param[2+2*nPrime]**2 *sigma(s, mPi, allowSubThr = allowSubThr)**3/sig**3
#	return retVal

INF = float('inf')
def dispersiveIntegral(s, param, renScale2, sCut, cutVal, allowSubThr = False):
	phiEval = phiV(s, param,  renScale2, allowSubThr = allowSubThr)
#	phiEvalMinDel = phiV(s, param,  renScale2, allowSubThr = allowSubThr)
#	phiEvalMacDel = phiV(s, param,  renScale2, allowSubThr = allowSubThr)
	mPi2    = mPi**2
	retVal  =- phiEval * constantRangeIntegral(s, 4*mPi2)
	retVal +=  constantRangeIntegral(s, sCut)*(phiEval - cutVal) # UpperLimit of \int_sThresh^sCut minus LowerLimit of \int_sCut^inf pi
	def kernel(sPrime):
		return (phiV(sPrime, param, renScale2, allowSubThr = allowSubThr)-phiEval)/(sPrime**3*(sPrime - s))
	retVal += integrate.quad(kernel, 4*mPi2, sCut)[0]
	return retVal

def constantRangeIntegral(s, sLimit):
	logArg  = 1.-s/sLimit
	retVal = clog(logArg)
	retVal *= 2*sLimit**2
	retVal += s*(2*sLimit + s)
	retVal /= 2*sLimit**2*s**3
	return retVal

class vectorFormFactor(parametersTrackingParameterization):
	def __init__(self, parameters, mu = .77**2, sCut = 4., cutVal = pi, x0 = 5.):
		self.x0          = x0
		self.mu          = mu
		self.sCut        = sCut
		self.cutVal      = cutVal
		self.nParAll     = 12
		self.allowSubThr = True
		self.makeLoadMap(parameters)
		self.piSheet     = 1
		self.Ksheet      = 1

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]
#		print par
		retVals = np.zeros((len(ms)), dtype = complex)
		for i, m in enumerate(ms):
			s = m**2
			retVals[i] = Fv(s, par, renScale2 = self.mu, sCut = self.sCut, cutVal = self.cutVal, x0 = self.x0, allowSubThr = self.allowSubThr, piSheet = self.piSheet, Ksheet = self.Ksheet).conjugate()
		return retVals

def main():
	params = [0.36674e-01, 0.31230e-02, .83386, .19771, 1.4974, .78518, 1.6855, .80109, 0.17254,-0.97591, 0.23374, 2.1982]
	name   = ["p_"+str(i) for i in range(len(params))]
	import ROOT

	nBins = 1000

	hist1 = ROOT.TH1D("h1", "h1", nBins, 0., 3.25)
	hist2 = ROOT.TH1D("h2", "h2", nBins, 0., 3.25)

	for i in range(nBins):
		s = hist1.GetXaxis().GetBinCenter(i+1)

		val1 = Fv(s, params, renScale2 = .77**2, sCut = 4, cutVal = pi, x0 = 5., allowSubThr = False)
		hist1.SetBinContent(i+1, abs(val1)**2)

		val2 = Fv(s, params, renScale2 = .77**2, sCut = 4, cutVal = pi, x0 = 5., allowSubThr = True)
		hist2.SetBinContent(i+1, abs(val2)**2)

#		print i, abs(val2)**2

#		print s,val1, val2
	hist1.SetMinimum(hist2.GetMinimum())
	hist1.Draw()
	hist2.Draw("SAME")
	raw_input()

if __name__ == "__main__":
	main()
