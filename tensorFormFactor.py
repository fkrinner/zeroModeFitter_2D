#!/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Python_ultra/Python-2.7.10/bin/python
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

from sheets import doCircularPlot

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

mPi = 0.13957018

def constantRangeIntegral(s, sLimit):
	logArg  = 1.-s/sLimit
#	retVal  = log(abs(logArg))
#	if logArg < 0.:
#		retVal += 1.j*pi
	retVal = clog(logArg)
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
		self.mPi         = mPi
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

	def FvTilde(s, param):
		part1 = (param[2]**2 + (param[8]*exp(1.j*param[9]) + param[10]*exp(1.j*param[11]))*s)/self.Den(s, param, 1)
		part2 = -param[ 8]*exp(1.j*param[ 9])*s/self.Den(s, param, 1)
		part3 = -param[10]*exp(1.j*param[11])*s/self.Den(s, param, 2)

