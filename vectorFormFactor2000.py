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

from parameterTrackingParameterizations import parametersTrackingParameterization, parameter
import scipy.integrate as integrate
import numpy as np

from sheets import doCircularPlot, findPoles
from plotVFFresults import parseFile, getNames

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

def constantRangeIntegral(s, sLimit):
	logArg  = 1.-s/sLimit
#	retVal  = log(abs(logArg))
#	if logArg < 0.:
#		retVal += 1.j*pi
	retVal  = clog(logArg)
	retVal *= 2*sLimit**2
	retVal += s*(2*sLimit + s)
	retVal /= 2*sLimit**2*s**3
	return retVal

def conformal(x,x0):
	return x0*atan(x/x0)

class vectorFormFactor(parametersTrackingParameterization):
	def __init__(self, parameters, mu = .77, sCut = mTau**2, cutVal = pi, x0 = 5.):
		self.x0          = x0
		self.mu          = mu
		self.sCut        = sCut
		self.cutVal      = cutVal
		self.nParAll     = 12
		self.allowSubThr = True
		self.piSheet     = 1
		self.Ksheet      = 1
		self.mPi         = mPi
		self.mK          = mK
		self.complexIntegration = False
		self.makeLoadMap(parameters)

	def __call__(self, ms, externalKinematicVariables = []):
		par = [p.value for p in self.parameters]
		retVals = np.zeros((len(ms)), dtype = complex)
		for i, m in enumerate(ms):
			s = m**2
			retVals[i] = self.Fv(s, par).conjugate()
		return retVals

	def sCall(self, s):
		par = [p.value for p in self.parameters]
		return self.Fv(s, par)

	def tildeCall(self, s):
		par = [p.value for p in self.parameters]
		return self.FvTilde(s, par)

	def Fv(self, s, param):
		mPi2     = self.mPi**2
		exponent =  param[0]*s/mPi2 + (param[1] - param[0]**2)*(s/mPi2)**2/2 + s**3/pi * self.dispersiveIntegral(s, param)
		if self.x0 is None:
			return exp(exponent)
		else:
			return exp(conformal(exponent.real, self.x0) + 1.j*exponent.imag)

	def dispersiveIntegral(self, s, param):
		phiEval = self.phiV(s, param)
		mPi2    = self.mPi**2
		retVal  =- phiEval * constantRangeIntegral(s, 4*mPi2)
		retVal +=  constantRangeIntegral(s, self.sCut)*(phiEval - self.cutVal) # UpperLimit of \int_sThresh^sCut minus LowerLimit of \int_sCut^inf pi
		if self.complexIntegration:
			def kernelReal(sPrime):
				return (self.phiV(sPrime, param)-phiEval)/(sPrime**3*(sPrime - s)).real
			retVal += integrate.quad(kernelReal, 4*mPi2, self.sCut)[0]
			def kernelImag(sPrime):
				return (self.phiV(sPrime, param)-phiEval)/(sPrime**3*(sPrime - s)).imag
			retVal += 1.j*integrate.quad(kernelImag, 4*mPi2, self.sCut)[0]

		else:
			def kernel(sPrime):
				return (self.phiV(sPrime, param)-phiEval)/(sPrime**3*(sPrime - s))
			retVal += integrate.quad(kernel, 4*mPi2, self.sCut)[0]
		return retVal

	def phiV(self, s,param):
		ampl = self.FvTilde(s, param)
		phi = atan(ampl.imag/ampl.real)
		if ampl.real < 0. or ampl.imag < 0.:
			phi += pi
		return phi

	def FvTilde(self, s, param):
		part1 = (param[2]**2 + (param[8]*exp(1.j*param[9]) + param[10]*exp(1.j*param[11]))*s)/(param[2]**2 - s+self.kappaRho(s,param) * (self.A(s, self.mPi, self.piSheet) + self.A(s, self.mK, self.Ksheet)/2))
		part2 = -param[ 8]*exp(1.j*param[ 9])*s/self.Den(s, param, 1)
		part3 = -param[10]*exp(1.j*param[11])*s/self.Den(s, param, 2)
		return part1 + part2 + part3

	def Den(self, s, param, nPrime):
		retVal = param[2+2*nPrime]**2 - s + self.kappaPrime(s, param, nPrime) * self.A(s, self.mPi, self.piSheet)
		return retVal

	def kappaRho(self, s,param):
		sig = self.sigma(param[2]**2, self.mPi)**3 + self.sigma(param[2]**2, self.mK)**3/2
		if sig == 0.:
			return 0.
		retVal = param[3]/param[2]*s/(pi*sig)
		return retVal

	def kappaPrime(self, s,param,nPrime):
		sig = self.sigma(param[2+2*nPrime]**2, self.mPi)
		if sig == 0.:
			return 0.
		retVal = param[3+2*nPrime]/param[2+2*nPrime] *s/(pi*sig**3)
		return retVal

	def sigma(self, s, mPart):
		xpr = 1. - 4*mPart**2/s + 0.j
		if xpr.real < 0. and not self.allowSubThr:
			return 0.
		return xpr**.5

	def A(self, s, mPart, sheet):
		sig    = self.sigma(s, mPart)
		retVal = log(mPart**2/self.mu**2) + 8*mPart**2/s - 5./3. + sig**3*clog((sig+1)/(sig-1.)) # +1.j pi to obtain the real part (at least above threshold)
		if sheet == 2:
			retVal += 2*pi* (4*mPart**2/s-1.+0.j)**1.5
		return retVal

def main():
	default_params =[parameter( .36674e-01, "lambP" ),
	                 parameter( .31230e-02, "lambPP"),
	                 parameter( .83386    , "M"     ),
	                 parameter( .19771    , "G"     ),
	                 parameter(1.4974     , "MP"    ),
	                 parameter( .78518    , "GP"    ),
	                 parameter(1.6855     , "MPP"   ),
	                 parameter( .80109    , "GPP"   ),
	                 parameter( .17254    , "alP"   ),
	                 parameter(-.97591    , "phP"   ),
	                 parameter( .23374    , "alPP"  ),
	                 parameter(2.1982     , "phPP"  )]

#	inFolder = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/global_vff_fits_subThresh/"
	inFolder = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/global_vff_fits/"
	names = [p.name for p in default_params]
	import ROOT
	hist = ROOT.TH2D("poles","poles", 100, 0.6, 1.6 ,100, 0., 1.)
	for fn in os.listdir(inFolder):
		if not '0123' in fn:
			continue
		if names is None:
			names = getNames(inFolder + fn)
		params = parseFile(inFolder + fn)
		pars = [parameter(params[i], names[i]) for i in range(12)]
		pars = default_params
		vff1 = vectorFormFactor(pars, sCut = 4)
		vff1.piSheet = 2
		vff1.Ksheet  = 1
		poles = findPoles(vff1.tildeCall, -.5, 2.5, -2., 2.)

		print pars[2].value**2,',',pars[2].value*pars[3].value
		print pars[4].value**2,',',pars[4].value*pars[5].value
		print pars[6].value**2,',',pars[6].value*pars[7].value

		for p in poles:
			print p.real**.5, p.imag/p.real**.5,p
			hist.Fill(p.real**.5, p.imag/p.real**.5)
		break
	hist.Draw("COLZ")
	raw_input()

if __name__ == "__main__":
	main()
