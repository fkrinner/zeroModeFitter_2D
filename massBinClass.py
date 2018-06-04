import pyRootPwa
import numpy as np
import numpy.linalg as la
import utils
from utils import getZeroHistSectors, normVector, getZeroModeNumber, printZeroStructure
import cmath
from cmath import phase, pi
from math import sin, cos
from modes import INTENS, PHASE, REAL, IMAG, INTENSNORM, INTENSTHEO, REALTHEO, IMAGTHEO, PHASETHEO, REIMCORRELATION
import sys

import scipy

def CMwrite(val):
	with open("comaLog", 'a') as outFile:
		outFile.write(str(val) + '\n')

class massBin:
	def __init__(self, nBin, realHists, imagHists, normHists, indexHists, coma, integralReal = None, integralImag = None):
		"""
		Initializer that sets data arrays
		"""
		self.bin3pi    = nBin
		self.binCenter = 0.52 + .04*nBin
		if not len(realHists) == len(imagHists) or not len(imagHists) == len(normHists) or not len(normHists) == len(indexHists):
			print "Numbers of histogams do not match:"
			print "  real:",len(realHists)
			print "  imag:",len(imagHists)
			print "  norm:",len(normHists)
			print " index:",len(indexHists)
			raise ValueError("Histogram size mismatch")
		self.nSect = len(realHists)
		if self.nSect == 0:
			raise ValueError("No histograms given.")
		self.nBins      = [ ]
		self.totalBins  =  0
		self.sectors    = [ ]
		for s in range(self.nSect):
			binMax = 0
			for bin in range(realHists[s].GetNbinsY()):
				m2Pi = realHists[s].GetYaxis().GetBinCenter( bin+1)
				m3Pi = realHists[s].GetXaxis().GetBinCenter(nBin+1)
				if utils.isValidPhaseSpace(m3Pi, m2Pi):
#				if realHists[s].GetBinContent(nBin + 1, bin+1) != 0.:
					binMax = bin
			self.nBins.append(binMax+1)
			self.totalBins += binMax+1
			self.sectors.append(realHists[s].GetTitle().split('_')[0])
		self.reals      = np.zeros((self.totalBins))
		self.imags      = np.zeros((self.totalBins))
		self.norms      = np.zeros((self.totalBins))
#	#	CMwrite("__init__")
		self.coma       = np.zeros((2*self.totalBins,2*self.totalBins))
		self.hasIntegralMatrix = False
		if integralReal and integralImag:
			self.hasIntegralMatrix = True
			self.integralMatrix    = np.zeros((self.totalBins, self.totalBins), dtype = complex)
		elif integralReal:
			raise RuntimeError("Cannot handle real integral matrix only, need also imaginary")
		elif integralImag:
			raise RuntimeError("Cannot handle imaginary integral matrix only, need also real")
		self.binCenters = np.zeros((self.totalBins))
		self.numLim     = 2.e-8
		self.ownPinv    = True
		count = 0
		for s in range(self.nSect):
			for bin in range(self.nBins[s]):
				self.reals[count]      = realHists[s].GetBinContent(nBin + 1, bin + 1)
				self.imags[count]      = imagHists[s].GetBinContent(nBin + 1, bin + 1)
				self.norms[count]      = normHists[s].GetBinContent(nBin + 1, bin + 1)
				self.binCenters[count] = realHists[s].GetYaxis().GetBinCenter(bin + 1)
				comaIndex = int(round(indexHists[s].GetBinContent(nBin + 1, bin + 1)))
				count2 = 0
				for s2 in range(self.nSect):
					for bin2 in range(self.nBins[s2]):
						comaIndex2 = int(round(indexHists[s2].GetBinContent(nBin + 1, bin2 + 1)))
						self.coma[2*count  , 2*count2  ] = coma.GetBinContent(2*comaIndex+1, 2*comaIndex2+1)
						self.coma[2*count  , 2*count2+1] = coma.GetBinContent(2*comaIndex+1, 2*comaIndex2+2)
						self.coma[2*count+1, 2*count2  ] = coma.GetBinContent(2*comaIndex+2, 2*comaIndex2+1)
						self.coma[2*count+1, 2*count2+1] = coma.GetBinContent(2*comaIndex+2, 2*comaIndex2+2)
						if self.hasIntegralMatrix:
							val = integralReal.GetBinContent(comaIndex+1, comaIndex2+1) + 1.j*integralImag.GetBinContent(comaIndex+1, comaIndex2+1)
							self.integralMatrix[count,count2] = val
						count2 += 1
				count +=1
		self.hasMassRange             = False
		self.makeComaInv()
		self.borders = [0]
		for i in range(self.nSect):
			self.borders.append(self.borders[-1] + self.nBins[i])
		self.nZero                    =  0
		self.zeroModes                = [ ]
		self.zeroModeNumbers          = [ ]
		self.zeroModeTitles           = [ ]
		self.zeroEigenvalues          = [ ]
		self.hasTheo                  = False
		self.chi2init                 = False
		self.zeroModesRemovedFromComa = False
		self.globalPhaseRemoved       = False
		self.specialCOMAs             = { }
		self.hasZMP                   = False
		self.zeroModeParameters       = None
		self.hasRandomizedAmplitudes  = False

	def getIntegralForFunction(self, sector, function):
		if not self.hasIntegralMatrix:
			raise RuntimeError("Cannot construct function integral without integral matrix")
		s        = self.getSectorInt(sector)
		startBin = self.borders[s  ]
		stopBin  = self.borders[s+1]
		ms       = self.binCenters[startBin:stopBin]
		vals     = function(ms, externalKinematicVariables = [self.binCenter])
		integral = 0.
		for i,b in enumerate(range(startBin,stopBin)):
			vals[i] *= self.norms[b]**.5
		for i,b in enumerate(range(startBin,stopBin)):
			for j,d in enumerate(range(startBin, stopBin)):
				integral += vals[i].conjugate()*vals[j]*self.integralMatrix[b,d]
		if integral.imag > self.numLim:
			raise ValueError("Integral witn non-zero imaginariy part encoutnered")
		return integral.real

	def unrandomize(self):
		"""
		Removes previously done gaussian randomization
		"""
		if not self.hasRandomizedAmplitudes:
			print "Warning: 'unrandomize' called, but has not been randomized"
		for i in range(self.totalBins):
			self.reals[i] -= self.randoms[2*i  ]
			self.imags[i] -= self.randoms[2*i+1]
		self.hasRandomizedAmplitudes = False
		del self.randoms

	def randomize(self):
		"""
		Randomize data points with a gaussian according to the covariance matrix
		"""
		if self.hasRandomizedAmplitudes:
			self.unrandomize()
		self.randoms = np.random.multivariate_normal(np.zeros((2*self.totalBins)), self.coma)
		for i in range(self.totalBins):
			self.reals[i] += self.randoms[2*i  ]
			self.imags[i] += self.randoms[2*i+1]
		self.hasRandomizedAmplitudes = True

	def getTheoryTotalMatrices(self, binRange = None):
		if not self.hasIntegralMatrix:
			raise RuntimeError("Cannot construct theory totals without integral matrix")
		totalMatrices = []
		for s in range(self.nSect):
			total = 0.
			startBin = self.borders[s  ]
			stopBin  = self.borders[s+1]
			masses = self.binCenters[startBin:stopBin]
			nFunc  = self.nFuncs[s]
			ampls  = [f(masses, externalKinematicVariables = [self.binCenter]) for f in self.funcs[s]]
			if binRange is not None:
				sect = self.sectors[s]
				if not sect in binRange:
					continue
				rang = range(binRange[sect][0], binRange[sect][1])
				for b in range(stopBin-startBin):
					if not b in rang:
						for f in range(len(ampls)):
							ampls[f][b] = 0.
			totalMatrix = np.zeros((nFunc,nFunc), dtype = complex)
			for iSect, iFull in enumerate(range(startBin,stopBin)):
				for jSect, jFull in enumerate(range(startBin,stopBin)):
					for iF in range(nFunc):
						for jF in range(nFunc):
							totalMatrix[iF,jF] += ampls[iF][iSect].conjugate() * self.integralMatrix[iFull,jFull] * ampls[jF][jSect]*(self.norms[iFull]*self.norms[jFull])**.5
			totalMatrices.append(totalMatrix)
		return totalMatrices

	def getSectorTotals(self, parameters, binRange = None, zeroModeSubtract = False):
		if not self.hasIntegralMatrix:
			raise RuntimeError("Cannot construct totals without integral matrix")
		ampls        = self.getCorrectedAmplitudes(parameters)
		if binRange is not None:
			for s in range(self.nSect):
				sect = self.sectors[s]
				if not sect in binRange:
					continue
				startBin = self.borders[s]
				stopBin  = self.borders[s+1]
				count    = 0
				rang     = range(binRange[sect][0], binRange[sect][1])
				for b in range(startBin, stopBin):
					if not count in rang:
						ampls[2*b  ] = 0.
						ampls[2*b+1] = 0.
					count += 1
		amplsComplex = np.zeros((len(ampls)/2), dtype = complex)
		realizedMatrix = np.zeros((len(ampls), len(ampls)))
		for i in range(len(ampls)/2):
			amplsComplex[i] = ampls[2*i] + 1.j*ampls[2*i+1]
			for j in range(len(ampls)/2):
				matrixValue =  self.integralMatrix[i,j]
				realizedMatrix[2*i  ,2*j  ] = matrixValue.real
				realizedMatrix[2*i  ,2*j+1] = matrixValue.imag
				realizedMatrix[2*i+1,2*j  ] =-matrixValue.imag
				realizedMatrix[2*i+1,2*j+1] = matrixValue.real
		realizedMatrix = realizedMatrix + np.transpose(realizedMatrix)
		totals = []
		for s in range(self.nSect):
			startBin = self.borders[s  ]
			stopBin  = self.borders[s+1]
			total = 0.
			amplsForErrors = np.zeros((len(ampls)))
			for i in range(startBin, stopBin):
				amplsForErrors[2*i  ] = ampls[2*i  ]
				amplsForErrors[2*i+1] = ampls[2*i+1]
				for j in range(startBin, stopBin):
					total += amplsComplex[i].conjugate() * self.integralMatrix[i,j] * amplsComplex[j]
				transformationMatrix = np.identity(len(self.coma))
				if zeroModeSubtract:
					transformationMatrix = np.zeros((len(self.coma),len(self.coma)))
					for z in range(self.nZero):
						z2 =0.
						for i in range(startBin, stopBin):
							if ampls[2*i] == 0.and ampls[2*i+1] == 0.:
								continue
							z2 += self.zeroModes[z][i]
						if z2 == 0.:
							continue
						for i in range(startBin, stopBin):
							if ampls[2*i] == 0. and ampls[2*i+1] == 0.:
								continue
							transformationMatrix[2*i  ,2*i  ] += 1.
							transformationMatrix[2*i+1,2*i+1] += 1.
							for j in range(startBin,stopBin):
								if ampls[2*j] == 0. and ampls[2*j+1] == 0.:
									continue
								transformationMatrix[2*i  ,2*j  ] -= self.zeroModes[z][i]*self.zeroModes[z][j]/z2
								transformationMatrix[2*i+1,2*j+1] -= self.zeroModes[z][i]*self.zeroModes[z][j]/z2
			if abs(total.imag) > self.numLim:
				raise ValueError("Total has non-vanishing imag part: " + str(total))
			jacobian = np.dot(realizedMatrix, amplsForErrors)
#			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#			print self.coma
#			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			err = np.dot(jacobian, np.dot(np.dot(np.dot(transformationMatrix,self.coma),transformationMatrix), jacobian))**.5
			totals.append((total, err))
#		print totals,'************************+*****************'
		return totals

	def initChi2(self, sectorFuncMap):
		"""
		Initializes the chi2 function
		"""
#		print "Set nPar to zero"
		self.nPar   =  0
		self.nFuncs = [ ]
		self.funcs  = self.toIntegerSectorMap(sectorFuncMap)
		for s in range(self.nSect):
			if not s in self.funcs:
				self.nFuncs.append(0)
			else:
				self.nFuncs.append(len(self.funcs[s]))
				for f in self.funcs[s]:
					self.nPar += f.nPar
#					print "Add",f.nPar,'to self.nPar:',self.nPar
		self.nFunc = 0
		for val in self.nFuncs:
			self.nFunc += val
		self.chi2init = True

	def chi2(self, pars = [], returnParameters = False):
		"""
		Evaluates the chi2 function, with shape parameters. If none are set, it uses, the internally set ones. non-shape parameters are calculated
		"""
		if not self.chi2init:
			raise RuntimeError("chi2 not inited, cannot evaluate")
			print pars
		if len(pars) is not 0:
			self.setShapeParameters(pars)
		A,B,C = self.getOwnTheoryABC()

#		print B
#		print A

		if self.ownPinv:
			zmPars  = -np.dot(B, utils.pinv(A + np.transpose(A), numLim = self.numLim))
		else:
			zmPars  = -np.dot(B, la.pinv(A + np.transpose(A)))


		chi2  = np.dot(zmPars,np.dot(A,zmPars)) + np.dot(zmPars,B) + C

		if returnParameters:
			return chi2, zmPars
		else:
			return chi2

	def fixedZMPchi2(self, pars):
		"""
		Returns a chi2 for the shape parameters and self.zeroModeParameters. The couplings are calculated.
		"""
		if not self.hasZMP and self.nZero > 0:
			raise RuntimeError("No zero mode parameters set")
		if pars is not None:
			self.setShapeParameters(pars)
		a,b,c = self.getOwnTheoryABC()
		A     = np.zeros((2*self.nFunc, 2*self.nFunc))
		B     = np.zeros((2*self.nFunc))
		C     = c
		for i in range(2*self.nZero):
			C    += b[i]*self.zeroModeParameters[i]
			for j in range(2*self.nZero):
				C += self.zeroModeParameters[i]*self.zeroModeParameters[j]*a[i,j]
		for i in range(2*self.nFunc):
			B[i] += b[2*self.nZero+i]
			for j in range(2*self.nZero):
				B[i] += (a[2*self.nZero+i,j]+a[j,2*self.nZero+i])*self.zeroModeParameters[j]
			for j in range(2*self.nFunc):
				A[i,j] += a[2*self.nZero + i, 2*self.nZero+j]
		if self.ownPinv:
			couplings = -np.dot(B, utils.pinv(np.transpose(A) + A, numLim = self.numLim))
		else:
			couplings = -np.dot(B, la.pinv(np.transpose(A) + A))
		return np.dot(couplings, np.dot(A,couplings)) + np.dot(B,couplings) + C

	def fixedZMPchi2_realCouplings(self, pars, phase, scan = True):
		"""
		Returns a chi2 for the shape parameters and self.zeroModeParameters. The real magnitudes of the couplings are calculated, they share a complex phase, which is given.
		"""
		if not self.hasZMP and self.nZero > 0:
			raise RuntimeError("No zero mode parameters set")
		if pars is not None:
			self.setShapeParameters(pars)
		
		a,b,c = self.getOwnTheoryABC()

		cosP  = cos(phase)
		sinP  = sin(phase)

		A     = np.zeros((self.nFunc, self.nFunc))
		B     = np.zeros((self.nFunc))
		C     = c
		for i in range(2*self.nZero):
			C    += b[i]*self.zeroModeParameters[i]
			for j in range(2*self.nZero):
				C += self.zeroModeParameters[i]*self.zeroModeParameters[j]*a[i,j]
		for i in range(self.nFunc):
			B[i] += cosP*b[2*self.nZero+2*i  ] 
			B[i] += sinP*b[2*self.nZero+2*i+1]
			for j in range(2*self.nZero):
				B[i] += cosP*(a[2*self.nZero+2*i  ,j]+a[j,2*self.nZero+2*i  ])*self.zeroModeParameters[j]
				B[i] += sinP*(a[2*self.nZero+2*i+1,j]+a[j,2*self.nZero+2*i+1])*self.zeroModeParameters[j]
			for j in range(self.nFunc):
				A[i,j] += cosP**2   * a[2*self.nZero + 2*i  , 2*self.nZero + 2*j  ]
				A[i,j] += cosP*sinP * a[2*self.nZero + 2*i  , 2*self.nZero + 2*j+1]
				A[i,j] += sinP*cosP * a[2*self.nZero + 2*i+1, 2*self.nZero + 2*j  ]
				A[i,j] += sinP**2   * a[2*self.nZero + 2*i+1, 2*self.nZero + 2*j+1]

		if self.ownPinv:
			couplings = -np.dot(B, utils.pinv(np.transpose(A) + A, numLim = self.numLim))
		else:
			couplings = -np.dot(B, la.pinv(np.transpose(A) + A))

		retVal = np.dot(couplings, np.dot(A,couplings)) + np.dot(B,couplings) + C

		return retVal

	def getFcnCplsABC(self, zeroModePars):
		a,b,c = self.getOwnTheoryABC()
		A     = np.zeros((2*self.nFunc, 2*self.nFunc))
		B     = np.zeros((2*self.nFunc))
		C     = c
		for i in range(2*self.nZero):
			C    += b[i]*zeroModePars[i]
			for j in range(2*self.nZero):
				C += zeroModePars[i]*zeroModePars[j]*a[i,j]
		for i in range(2*self.nFunc):
			B[i] += b[2*self.nZero+i]
			for j in range(2*self.nZero):
				B[i] += (a[2*self.nZero+i,j]+a[j,2*self.nZero+i])*zeroModePars[j]
			for j in range(2*self.nFunc):
				A[i,j] += a[2*self.nZero + i, 2*self.nZero+j]
		return A,B,C

	def getFcnCplsHess(self, zeroModePars):
		A,B,C = self.getFcnCplsABC(zeroModePars)
		return np.transpose(A) + A

	def getFcnCpls(self, zeroModePars):
		A,B,C = self.getFcnCplsABC(zeroModePars)
		if self.ownPinv:
			couplings = -np.dot(B, utils.pinv(np.transpose(A) + A, numLim = self.numLim))
		else:
			couplings = -np.dot(B, la.pinv(np.transpose(A) + A))
		return couplings

	def getNonShapeUncertainties(self, pars = []):
		"""
		Gets the unceratinties on non-shape parameters at the minimum of chi2(...)
		"""
		if not self.chi2init:
			raise RuntimeError("Chi2 not initialized, cannot evaluate")
		if len(pars) > 0:
			self.setShapeParameters(pars)
		A,B,C = self.getOwnTheoryABC()
		unc = []
		for i in range(len(A)):
			unc.append(A[i,i]**-.5)
		return unc

	def getNonShapeParameters(self, pars = []):
		"""
		Gets the non-shape parameters at the minimum of chi2(...)
		"""
		if not self.chi2init:
			raise RuntimeError("Chi2 not initialized, cannot evaluate")
		if len(pars) > 0:
			self.setShapeParameters(pars)
		A,B,C = self.getOwnTheoryABC()
		if self.ownPinv:
			return  -np.dot(B, utils.pinv(A + np.transpose(A), numLim = self.numLim))
		else:
			return  -np.dot(B, la.pinv(A + np.transpose(A)))

	def getChi2forNonShapeParameters(self,nonShapeParameters, shapePars = []):
		"""
		Gets the chi2 for a given set of non-shape parameters
		"""
		if not self.chi2init:
			raise RuntimeError("Chi2 not initialized, cannot evaluate")
		if len(shapePars) > 0:
			self.setShapeParameters(shapePars)
		A,B,C = self.getOwnTheoryABC()
		return np.dot(nonShapeParameters, np.dot(A,nonShapeParameters)) + np.dot(B,nonShapeParameters) + C

	def compareTwoZeroModeCorrections(self, params1, params2):
		"""
		Compares two different corrections of the zero-modes
		"""
		if self.zeroModesRemovedFromComa:
			raise RuntimeError("Comparison does not make sense with removed zero-mode directions")
		if not len(params1) == 2*self.nZero:
			raise ValueError("Number of zero modes does not match (1)")
		if not len(params2) == 2*self.nZero:
			raise ValueError("Number of zero modes does not match (2)")
		deltas = np.zeros((2*self.totalBins))
		print '--------in mass bin class'
		print '--------',params1
		print '--------',params2
		print '--------out mass bin class'
		for b in range(self.totalBins):
			for z in range(self.nZero):
				zm = self.zeroModes[z][b]
#				coeff = params1[z] - params2[z]
				deltas[2*b  ] = zm*(params1[2*z  ] - params2[2*z  ])
				deltas[2*b+1] = zm*(params1[2*z+1] - params2[2*z+1])
		print '-------->>>',np.dot(deltas, deltas)
		return np.dot(deltas, np.dot(self.comaInv, deltas))

	def nParAll(self):
		"""
		Total number of paramters
		"""
		return 2*self.nZero + 2*self.nFunc + self.nPar

	def phaseChi2(self, pars):
		"""
		Gets the chi2 for the PHASE method
		"""
		return self.modeChi2(pars, PHASE)

	def intensChi2(self, pars):
		"""
		Gets the chi2 for the INTENS method
		"""
		return self.modeChi2(pars, INTENS)

	def modeChi2(self, pars, mode = PHASE):
		"""
		Gets the chi2 for a special mode, PHASE of INTENS
		"""
		if not self.chi2init:
			raise RuntimeError("chi2 not inited, cannot evaluate")
		if not len(pars) == self.nParAll():
			raise ValueError("Number of parameters does not match "+ str(len(pars)) + " vs " + str(self.nParAll()))
		corrAmpl = self.getCorrectedAmplitudes(pars[:2*self.nZero])
		self.setShapeParameters(pars[2*(self.nZero + self.nFunc):])
		self.setTheoryFromOwnFunctions(pars[2*self.nZero:2*(self.nZero + self.nFunc)])

		deltas = np.zeros((self.totalBins))
		for s in range(self.nSect):
			if not s in self.funcs:
				continue
			startBin = self.borders[s  ]
			stopBin  = self.borders[s+1]
			for bin in range(startBin, stopBin):
				if self.hasMassRange:
					if s in self.massRanges:
						binCenterMass = self.binCenters[bin]
						if binCenterMass < self.massRanges[s][0] or binCenterMass >= self.massRanges[s][1]:
							continue
				if mode == PHASE:
#					print phase(self.theo[bin])
					deltas[bin] = phase(corrAmpl[2*bin] + 1.j*corrAmpl[2*bin+1]) - phase(self.theo[bin])
					if deltas[bin] > pi:
						deltas[bin] -= 2*pi
					if deltas[bin] < -pi:
						deltas[bin] += 2*pi
#					print phase(corrAmpl[2*bin] + 1.j*corrAmpl[2*bin+1]),'-',phase(self.theo[bin]),'=',deltas[bin]

				elif mode == INTENS:
					deltas[bin] = corrAmpl[2*bin]**2 + corrAmpl[2*bin+1]**2 - abs(self.theo[bin])**2
		coma = self.getSpecialComaInv(mode)
#		printZeroStructure(coma)
		retVal = np.dot(deltas, np.dot(coma, deltas))
#		print np.dot(deltas, deltas),"{{{{{{{{{{{}"

		removePhaseDirection = True
		if removePhaseDirection:
			entry          = (1./len(deltas))**.5 # The direction of a global phase rotation is (1., 1., 1., 1., ...) / len(deltas)**.5
			sp             = 0.
			for i in range(len(deltas)):
				sp += entry * deltas[i]
			retVal += sp**2*retVal
		return retVal

	def getSpecialComaInv(self, mode = PHASE):
		"""
		Returns the transformed covarinace matrix for a special method. Stores it also, to do not have to do this over and over again
		"""
		if mode in self.specialCOMAs:
			return self.specialCOMAs[mode]
		mapy = np.copy(self.coma) # coMA co PY
		if self.hasMassRange:
			zeroIndices = []
			for s in range(self.nSect):
				for i,b in enumerate(range(self.borders[s],self.borders[s+1])):
					if s in self.massRanges:
						binCenterMass = self.binCenters[b]
						if binCenterMass < self.massRanges[s][0] or binCenterMass >= self.massRanges[s][1]:
							zeroIndices.append(b)
			for i in range(len(mapy)):
				for zi in zeroIndices:
					mapy[i,2*zi  ] = 0.
					mapy[i,2*zi+1] = 0.
					mapy[2*zi  ,i] = 0.
					mapy[2*zi+1,i] = 0.

		jacobian = np.zeros((self.totalBins, 2*self.totalBins))
		for i in range(self.totalBins):
			if mode == INTENS:
				jacobian[i, 2*i  ] = 2*self.reals[i]
				jacobian[i, 2*i+1] = 2*self.imags[i]
			elif mode == PHASE:
				re = self.reals[i]
				im = self.imags[i]
				if re == 0.:
					if im > 0.:
						jacobian[i,2*i  ] = -1./im
					else:
						jacobian[i,2*i  ] = 1./im
				else:
					common = 1. + im**2/re**2
					jacobian[i,2*i  ] = -im/re**2/common
					jacobian[i,2*i+1] = 1./re/common
		if self.onwPinv:
			retVal = utils.pinv(np.dot(jacobian, np.dot(mapy, np.transpose(jacobian))), numLim = self.numLim)
		else:
			retVal = la.pinv(np.dot(jacobian, np.dot(mapy, np.transpose(jacobian))))
		self.specialCOMAs[mode] = retVal
		return retVal

	def setShapeParameters(self, pars):
		"""
		Sets the shape parameters of the functions
		"""
		if not self.chi2init:
			raise RuntimeError("chi2 not inited, no knowledge about shape parameters")
		if not len(pars) == self.nPar:
			print "GRAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
			print pars
			print "GRAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
			raise ValueError("Number of shape parameters does not match len(pars) = " + str(len(pars)) + " != self.nPar = " + str(self.nPar))
		countPar = 0
		for s in range(self.nSect):
			if not s in self.funcs:
				continue
			for f in self.funcs[s]:
				f.setParameters(pars[countPar:countPar + f.nPar])
				countPar += f.nPar

	def eigenbasisChi2contributions(self, params):
		"""
		Gets the chi2-contributions in terms of eigenvectors of the covariance matrix
		"""
		count = 0
		deltas = self.getCorrectedAmplitudes(params[:2*self.nZero])
		self.setTheoryFromOwnFunctions(params[2*self.nZero:])
		for i in range(self.totalBins):
			deltas[2*i  ] -= self.theo[i].real
			deltas[2*i+1] -= self.theo[i].imag
		vals, vecs = la.eig(self.coma)
		retVal = []
		for i,val in enumerate(vals):
			proj = np.dot(vecs.T[i], deltas)
			retVal.append((val, proj))
		return retVal		

	def setZeroModeParameters(self, zmp):
		"""
		Sets fixed zero mode parameters
		"""
		if not len(zmp) == 2*self.nZero:
			raise IndexError("Number of zero mode parameters does not match")
		self.hasZMP             = True
		self.zeroModeParameters = zmp

	def setParametersAndErrors(self, pars, errors):
		"""
		Sets parameters and errors of the functions
		"""
		if not self.chi2init:
			print "chi2 not inited, no knowledge about shape parameters"
			return False
		if not len(errors) == self.nPar:
			print "Number of parameter errors does not match len(errors) = " + str(len(errors)) + " != self.nPar = " + str(self.nPar)
			return False
		if not len(pars) == self.nPar:
			print "Number of shape parameters does not match len(pars) = " + str(len(pars)) + " != self.nPar = " + str(self.nPar)
			return False
		countPar = 0
		for s in range(self.nSect):
			if not s in self.funcs:
				continue
			for f in self.funcs[s]:
				if not f.setParametersAndErrors(pars[countPar:countPar + f.nPar], errors[countPar:countPar + f.nPar]):
					print "Could not set parameter in function"
					return False
				countPar += f.nPar
		return True

	def getOwnTheoryABC(self):
		"""
		Gets the A B and C from the own theory functions: chi2(p) = pAP + Bp + C
		"""
		if not self.chi2init:
			raise RuntimeError("chi2 not inited, cannot get ABC")
		nTotal = self.nZero + self.nFunc
		A = np.zeros((2*nTotal, 2*nTotal))
		B = np.zeros((2*nTotal))
		C = 0.
		countNfunc = self.nZero
		for s in range(self.nSect):
			if not s in self.funcs:
				continue
			nFunc  = self.nFuncs[s]
			masses = self.binCenters[self.borders[s]:self.borders[s+1]]
			ampls  = [f(masses, externalKinematicVariables = [self.binCenter]) for f in self.funcs[s]]
			massRange = None
			if self.hasMassRange:
				if s in self.massRanges:
					massRange = self.massRanges[s]
			a,b,c  = self.getTheoryABC(s, ampls, massRange = massRange)
			for i in range(len(b)):
				if not b[i].imag == 0.:
					if abs(b[i].imag) > self.numLim:
						raise ValueError("Complex value in b " + str(b[i]))
					else:
						b[i] = b[i].real
				for j in range(len(b)):
					if not a[i,j].imag == 0.:
						if abs(a[i,j].imag) > self.numLim:
							raise ValueError("Complex value in a " + str(a[i,j]))
						else:
							a[i,j] = a[i,j].real
			for i in range(self.nZero):
				B[2*i  ] += b[2*i  ].real
				B[2*i+1] += b[2*i+1].real
				for j in range(self.nZero):
					A[2*i  ,2*j  ] += a[2*i  ,2*j  ].real
					A[2*i  ,2*j+1] += a[2*i  ,2*j+1].real
					A[2*i+1,2*j  ] += a[2*i+1,2*j  ].real
					A[2*i+1,2*j+1] += a[2*i+1,2*j+1].real
				for j in range(nFunc):
					A[2*i  ,2*(j+countNfunc)  ] += a[2*i  ,2*(self.nZero+j)  ].real
					A[2*i  ,2*(j+countNfunc)+1] += a[2*i  ,2*(self.nZero+j)+1].real
					A[2*i+1,2*(j+countNfunc)  ] += a[2*i+1,2*(self.nZero+j)  ].real
					A[2*i+1,2*(j+countNfunc)+1] += a[2*i+1,2*(self.nZero+j)+1].real

					A[2*(j+countNfunc)  ,2*i  ] += a[2*(self.nZero+j)  ,2*i  ].real
					A[2*(j+countNfunc)  ,2*i+1] += a[2*(self.nZero+j)  ,2*i+1].real
					A[2*(j+countNfunc)+1,2*i  ] += a[2*(self.nZero+j)+1,2*i  ].real
					A[2*(j+countNfunc)+1,2*i+1] += a[2*(self.nZero+j)+1,2*i+1].real
			for i in range(nFunc):
				B[2*(i+countNfunc)  ] += b[2*(self.nZero+i)  ].real
				B[2*(i+countNfunc)+1] += b[2*(self.nZero+i)+1].real
				for j in range(nFunc):
					A[2*(i+countNfunc)  ,2*(j+countNfunc)  ] += a[2*(self.nZero+i)  ,2*(self.nZero+j)  ].real
					A[2*(i+countNfunc)  ,2*(j+countNfunc)+1] += a[2*(self.nZero+i)  ,2*(self.nZero+j)+1].real
					A[2*(i+countNfunc)+1,2*(j+countNfunc)  ] += a[2*(self.nZero+i)+1,2*(self.nZero+j)  ].real
					A[2*(i+countNfunc)+1,2*(j+countNfunc)+1] += a[2*(self.nZero+i)+1,2*(self.nZero+j)+1].real
			C += c
			countNfunc += nFunc
		return A,B,C

	def setZeroTheory(self):
		"""
		Sets the own theory evaluation to zero
		"""
		self.theo = np.zeros((self.totalBins), dtype = complex)

	def setTheory(self, sectParFunctMap = {}):
		"""
		Sets the own theory evaluations to the given or set functions
		"""
		spfm = self.toIntegerSectorMap(sectParFunctMap)

		self.theo = np.zeros((self.totalBins), dtype = complex)
		if spfm == {}:
			print "setTheory called with no functions. Setting zero theory as dummy..."
		for s in range(self.nSect):
			if not s in spfm:
				continue
			startBin = self.borders[s  ]
			stopBin  = self.borders[s+1]
			binning  = self.binCenters[startBin:stopBin]
			ampl = np.zeros((len(binning)), dtype = complex)
			for pf in spfm[s]:
				if len(pf) == 2: # If shape paremeters of the functions are to be used
					ampl += pf[0] * pf[1](binning, externalKinematicVariables = [self.binCenter])
				elif len(pf) == 3: # If shape parameters are given as third element
					ampl += pf[0] * pf[1](binning, pf[2],externalKinematicVariables = [self.binCenter] )
				else:
					raise ValueError("Format of sectParFunctMap is wrong")
			for i,b in enumerate(range(startBin, stopBin)):
				self.theo[b] = ampl[i]*self.norms[b]**.5
				if self.hasMassRange:
					if s in self.massRanges:
						binCenterMass = self.binCenters[b]
						if binCenterMass < self.massRanges[s][0] or binCenterMass >= self.massRanges[s][1]:
							self.theo[b] = 0.+0.j
		self.hasTheo = True

	def unifyComa(self):
		"""
		Sets unity as covariance matrix
		"""
#	#	CMwrite("unifyComa")
		for i in range(len(self.coma)):
			for j in range(len(self.coma)):
				if i == j and not self.coma[i,i] == 0.:
					self.coma[i,j] = 1.
					self.comaInv[i,j] = 1.
				else:
					self.coma[i,j] = 0.
					self.comaInv[i,j] = 0.
		self.specialCOMAs = {}

	def makeComaInv(self):
		"""
		Inverts the covariance matrix (To be able to change the inversion method globally)
		"""
#	Change here: Check sector range map an invert submatrix
#
#
		mapy = np.copy(self.coma) # coMA co PY
		if self.hasMassRange:
			zeroIndices = []
			for s in range(self.nSect):
				for i,b in enumerate(range(self.borders[s],self.borders[s+1])):
					if s in self.massRanges:
						binCenterMass = self.binCenters[b]
						if binCenterMass < self.massRanges[s][0] or binCenterMass >= self.massRanges[s][1]:
							zeroIndices.append(b)
			for i in range(len(self.coma)):
				for zi in zeroIndices:
					mapy[i,2*zi  ] = 0.
					mapy[i,2*zi+1] = 0.
					mapy[2*zi  ,i] = 0.
					mapy[2*zi+1,i] = 0.
		if self.ownPinv:
			self.comaInv = utils.pinv(mapy, numLim = self.numLim)
		else:
			self.comaInv = la.pinv(mapy)
		

	def removeAllCorrelations(self, removeReImCorrel = True):
		"""
		Removes ALL correlations from the covariance matrix
		"""
		dim = len(self.coma)/2
#	#	CMwrite("removeAllCorrelations")
		for i in range(dim):
			for j in range(dim):
				if not i == j:
					self.coma[2*i  ,2*j  ] = 0.		
					self.coma[2*i+1,2*j  ] = 0.
					self.coma[2*i  ,2*j+1] = 0.
					self.coma[2*i+1,2*j+1] = 0.
				elif removeReImCorrel:
					self.coma[2*i+1,2*j  ] = 0.
					self.coma[2*i  ,2*j+1] = 0.
		self.makeComaInv()
		self.specialCOMAs = {}

	def removeZeroModeFromComa(self):
		"""
		Removes the dorection of the zero-mode from the covariance matrix
		"""
#		print 'Remove zm'
		if self.zeroModesRemovedFromComa:
			print "DO NOT REMOVE ZER MODE A SECOND TIME"
			return

		self.zeroModesRemovedFromComa = True
		if len(self.zeroModes) == 0:
			return
		dim = len(self.coma)
		transformationMatrix = np.identity(dim)
#	#	CMwrite("removeZeroModeFromComa")
		for z in range(self.nZero):
			for i in range(dim/2):
				for j in range(dim/2):
					transformationMatrix[2*i  ,2*j  ] -= self.zeroModes[z][i] * self.zeroModes[z][j]
					transformationMatrix[2*i+1,2*j+1] -= self.zeroModes[z][i] * self.zeroModes[z][j]
		self.coma         = np.dot(transformationMatrix, np.dot(self.coma, transformationMatrix)) # no transpose needed, since transformationMatrix is symmetric
		self.makeComaInv()
		self.specialCOMAs = {}

	def addComaValueForZeroMode(self, val, unitsOf = 'smallestComaValue'):
		if not self.zeroModesRemovedFromComa:
			raise RuntimeError("Call 'removeZeroModeFromComa()' first.")
		if unitsOf == 'smallestComaValue':
			vals, vecs = la.eig(self.coma)
			minVal = float('inf')
			for v in vals:
				if not v < self.numLim:
					minVal = min(v, minVal)
			value = val * minVal
		elif unitsOf == "one":
			value = val
		else:
			raise ValueError("Unknown option unitsOf = '" + untisOf + "'")
		dim = len(self.coma)
#	#	CMwrite("addComaValueForZeroMode("+str(val)+', '+unitsof+')')
		for z in range(self.nZero):
			for i in range(dim/2):
				for j in range(dim/2):
					self.coma[2*i  ,2*j  ] += value* self.zeroModes[z][i] * self.zeroModes[z][j]
					self.coma[2*i+1,2*j+1] += value* self.zeroModes[z][i] * self.zeroModes[z][j]
		self.makeComaInv()
		self.specialCOMAs = {}

	def zeroModeMultiplicationCheck(self, coeff = 1.):
		"""
		Multiplies the inverted coma with the zero-modes and prints the result
		"""
		dim = len(self.comaInv)
		for z in range(self.nZero):
			val = 0.
			for i in range(dim/2):
				for j in range(dim/2):
					val += self.zeroModes[z][i]*self.comaInv[2*i  ,2*j  ]*self.zeroModes[z][j]
					val += self.zeroModes[z][i]*self.comaInv[2*i  ,2*j+1]*self.zeroModes[z][j]
					val += self.zeroModes[z][i]*self.comaInv[2*i+1,2*j  ]*self.zeroModes[z][j]
					val += self.zeroModes[z][i]*self.comaInv[2*i+1,2*j+1]*self.zeroModes[z][j]
			print "Check for zero mode",str(z),":",str(val)


	def getGlobalPhaseDirection(self):
		"""
		Gets the direction of a global phase rotation
		"""
		retVal = np.zeros((2*self.totalBins))
		for i in range(self.totalBins):
			retVal[2*i  ] = -self.imags[i]
			retVal[2*i+1] =  self.reals[i]

		return normVector(retVal)

	def removeGlobalPhaseFromComa(self):
		"""
		Removes the direction of a global phase rotation from the covariance matrix
		"""
		print "Phase removed <---!"
		if self.globalPhaseRemoved:
			print "DO NOT REMOVE THE GLOBAL PHASE A SECOND TIME"
			return 
		self.globalPhaseRemoved = True
		phaseDirection = self.getGlobalPhaseDirection()
		dim = 2*self.totalBins
		transformationMatrix = np.identity(dim)
#	#	CMwrite("removeGlobalPhaseFromComa")
		for i in range(dim):
			for j in range(dim):
				transformationMatrix[i,j] -= phaseDirection[i] * phaseDirection[j]
		self.coma         = np.dot(transformationMatrix, np.dot(self.coma, transformationMatrix)) # no transpose needed, since transformationMatrix is symmetric
		self.makeComaInv()
		self.specialCOMAs = {}

	def setMassRanges(self, massRanges):
		"""
		Sets the mass ranges for the chi2 functions
		"""
		if self.hasMassRange:
			raise RuntimeError("Cannot setMassRange(...) twice, since the COMA has to be cut")
		self.massRanges   = self.toIntegerSectorMap(massRanges)
#		zeroIndices = []
#	#	CMwrite("setMassRanges("+str(massRanges)+')')
	#	for s in range(self.nSect):
	#		for i,b in enumerate(range(self.borders[s],self.borders[s+1])):
	#			if s in self.massRanges:
	#				binCenterMass = self.binCenters[b]
	#				if binCenterMass < self.massRanges[s][0] or binCenterMass >= self.massRanges[s][1]:
	#					zeroIndices.append(b)
	#	for i in range(len(self.coma)):
	#		for zi in zeroIndices:
	#			self.coma[i,2*zi  ] = 0.
	#			self.coma[i,2*zi+1] = 0.
	#			self.coma[2*zi  ,i] = 0.
	#			self.coma[2*zi+1,i] = 0.

		self.hasMassRange = True
		self.makeComaInv()
		self.specialCOMAs = {}

	def setTheoryFromOwnFunctions(self, parameterList, restrictToRange = True):
		"""
		Sets theory curves from own functions
		"""
		if self.hasMassRange and not restrictToRange:
			raise RuntimeError("Setting a mass range now changes the COMA, this combination of features is not avalable anymore")
		if not self.chi2init:
			raise RuntimeError("Chi2 not inited, does not have own functions")
		if not len(parameterList) == 2*self.nFunc:
			raise ValueError("Number of parameters does not match "+ str(len(parameterList)) + ' ' + str(2*self.nFunc))

		countFunc = 0
		self.theo = np.zeros((self.totalBins), dtype = complex)
		for s in range(self.nSect):
			if not s in self.funcs:
				continue
			nFunc = self.nFuncs[s]
			masses = self.binCenters[self.borders[s]:self.borders[s+1]]

			ampl = np.zeros((len(masses)), dtype = complex)
			for i,f in enumerate(self.funcs[s]):
				cpl = parameterList[2*(i+countFunc)] + 1.j*parameterList[2*(i+countFunc)+1]
				ampl += cpl * f(masses, externalKinematicVariables = [self.binCenter])
			countFunc += nFunc
			for i,b in enumerate(range(self.borders[s],self.borders[s+1])):
				self.theo[b] = ampl[i]*self.norms[b]**.5
				if self.hasMassRange and restrictToRange:
					if s in self.massRanges:
						binCenterMass = self.binCenters[b]
						if binCenterMass < self.massRanges[s][0] or binCenterMass >= self.massRanges[s][1]:
							self.theo[b] = 0.+0.j
		self.hasTheo = True

	def getNDF(self):
		"""
		Return NDF as (NDF, 2*nZero, 2*nFunc, nPar)
		"""
		NDF = 0
		for s in range(self.nSect):
			if not s in self.funcs:
				continue
			for i,b in enumerate(range(self.borders[s],self.borders[s+1])):
				if self.hasMassRange:
					if s in self.massRanges:
						binCenterMass = self.binCenters[b]
						if binCenterMass < self.massRanges[s][0] or binCenterMass >= self.massRanges[s][1]:
							continue
				NDF += 2
		return (NDF, 2*self.nZero, 2*self.nFunc, self.nPar)

	def writeZeroModeCoefficients(self, coefficients, outFileName, tBin = "<tBin>", mode = 'a'):
		"""
		Writes the zero-mode coefficients to a text file
		"""
		tBin  = str(tBin)
		cmplx = False
		if len(coefficients) == self.nZero:
			cmplx = True
		elif not len(coefficients) == 2*self.nZero:
			raise IndexError("Number of coefficients for "+str(self.nZero)+" modes does not match (neither real nor complex): " +str(len(coefficients)))
		with open(outFileName, mode) as out:
			out.write(tBin + ' ' + str(self.bin3pi))
			for i in range(self.nZero):
				out.write(' ' + self.zeroModeTitles[i])
				if cmplx:
					out.write(' ' + str(coefficients[i].real) + ' ' + str(coefficients[i].imag))
				else:
					out.write(' ' + str(coefficients[2*i]) + ' ' + str(coefficients[2*i+1]))
			out.write('\n')

	def renormZeroModes(self):
		"""
		Renors the zero modes
		"""
		for z in range(self.nZero):
			for i in range(self.totalBins):
				self.zeroModes[z][i] *= self.norms[i]**.5

	def addZeroMode(self, modeBorders, modeHist, eigenvalueHist = None):
		"""
		Add a zero mode, determines, whether it is needed and sets it accordingly
		"""
		newMode  = np.zeros((self.totalBins))
		modeList = getZeroHistSectors(modeHist)
		modeMap  = {}
		if not eigenvalueHist == None:
			if not len(self.zeroModes) == len(self.zeroEigenvalues):
				print "ERROR: Eigenvalue hist given, but sizes of values and vectors do not match"
				return False
			self.zeroEigenvalues.append(eigenvalueHist.GetBinContent(self.bin3pi + 1))
		else:
			if not len(self.zeroEigenvalues) == 0:
				print "ERROR: Eigenvalue was previously added, but is now missing"
				return False
		for s, sect in enumerate(self.sectors):
			for S, SECT in enumerate(modeList):
				if sect == SECT:
					modeMap[s] = S
		for s in range(self.nSect):
			if not s in modeMap:
				continue
			nMode = modeMap[s]
			count = 0
			for i in range(self.borders[s], self.borders[s+1]):
				newMode[i] = modeHist.GetBinContent(self.bin3pi + 1, modeBorders[nMode] + count + 1)
				count += 1
		allZero = True
		for i in range(self.totalBins):
			if not newMode[i] == 0.:
				allZero = False
				break
		if allZero:
#			print "WARNING: Zero mode is vanishing everywhere, do not add"
			return False
		self.zeroModeNumbers.append(getZeroModeNumber(modeHist))
		self.zeroModeTitles.append(modeHist.GetTitle())
		self.zeroModes.append(newMode)
		self.nZero = len(self.zeroModes)

		return True

	def getTheoryABC(self, sector, parametrizations, massRange = None):
		"""
		Gets the A B and C for evalueations from the set theory: chi2(p) = pAP + Bp + C
		"""
		nPara   = len(parametrizations)
		if sector < 0 or sector >= self.nSect:
			raise IndexError("Sector "+str(sector)+" invalid")
		nBinMin = self.borders[sector  ]
		nBinMax = self.borders[sector+1]
		nBins   = nBinMax - nBinMin
		for p in range(nPara):
			if not len(parametrizations[p]) == nBins:
				raise IndexError("Prameterization has the wrong size")
		Ampls = np.zeros((2*self.totalBins))
		ZP    = np.zeros((2*(self.nZero + nPara), 2*self.totalBins))
		count = 0
		for bin in range(nBinMin, nBinMax):
			if massRange:
				mass = self.binCenters[bin]
				if mass < massRange[0] or mass >= massRange[1]:
					count += 1
					continue
			Ampls[2*bin  ] = self.reals[bin]
			Ampls[2*bin+1] = self.imags[bin]
			for z in range(self.nZero):
				ZP[2*z  ,2*bin  ] = self.zeroModes[z][bin]
				ZP[2*z+1,2*bin+1] = self.zeroModes[z][bin]
			for p in range(nPara):
				ampl = parametrizations[p][count] * self.norms[bin]**.5
	#			print ampl
				ZP[2*(p+self.nZero)  ,2*bin  ] = - ampl.real # Minus sign from A + C_z Z - C_p P, plus sign from    Re = R*r - I*i
				ZP[2*(p+self.nZero)+1,2*bin  ] =   ampl.imag # Minus sign from A + C_z Z - C_p P, second minus from Re = R*r - I*i
				ZP[2*(p+self.nZero)  ,2*bin+1] = - ampl.imag # Minus sign from A + C_z Z - C_p P, plus sign from    Im = R*i + I*r
				ZP[2*(p+self.nZero)+1,2*bin+1] = - ampl.real # Minus sign from A + C_z Z - C_p P, plus sign from    Im = R*i + I*r
			count += 1
		AC = np.dot(Ampls, self.comaInv)
		C  = np.dot(Ampls, AC)
		B  = 2*np.dot(ZP, AC)
		A  = np.dot(ZP, np.transpose(np.dot(ZP, self.comaInv)))
		return A,B,C

	def getSmoothnessABC(self,other):
		"""
		Gets the A B and C for the smoothness method
		"""
		if not self.sectors == other.sectors:
			raise ValueError("Sectors do not match")
		totalNzero = self.nZero + other.nZero

		Dself  = np.zeros((2*self.totalBins))
		Dother = np.zeros((2*other.totalBins))
		Zself  = np.zeros((2*totalNzero, 2*self.totalBins))
		Zother = np.zeros((2*totalNzero, 2*other.totalBins))
		for s in range(self.nSect):
			startSelf  = self.borders[s]
			startOther = other.borders[s]
			nBins      = min(self.borders[s+1] - self.borders[s], other.borders[s+1] - other.borders[s])
			for i in range(nBins):
				delRe = self.reals[startSelf+i] - other.reals[startOther+i]
				delIm = self.imags[startSelf+i] - other.imags[startOther+i]
				Dself[2*(startSelf + i)  ] = delRe
				Dself[2*(startSelf + i)+1] = delIm
				Dother[2*(startOther + i)  ] = delRe
				Dother[2*(startOther + i)+1] = delIm
				for z in range(self.nZero):
					zeroVal = self.zeroModes[z][startSelf+i]
					Zself[2*z  ,2*(startSelf + i)    ] = zeroVal
					Zself[2*z+1,2*(startSelf + i)+1  ] = zeroVal
					Zother[2*z  ,2*(startOther + i)  ] = zeroVal
					Zother[2*z+1,2*(startOther + i)+1] = zeroVal
				z0 = self.nZero
				for z in range(other.nZero):
					zeroVal = -other.zeroModes[z][startOther+i]
					Zself[2*(z+z0)   , 2*(startSelf + i)  ] = zeroVal
					Zself[2*(z+z0)+1 , 2*(startSelf + i)+1] = zeroVal
					Zother[2*(z+z0)  , 2*(startOther+ i)  ] = zeroVal
					Zother[2*(z+z0)+1, 2*(startOther+ i)+1] = zeroVal
		DCself  = np.dot(Dself,  self.comaInv)
		DCother = np.dot(Dother, other.comaInv)

		Cself  = np.dot(Dself,  DCself)
		Cother = np.dot(Dother, DCother)
		Bself  = np.dot(Zself,  DCself)
		Bother = np.dot(Zother, DCother)
		Aself  = np.dot(Zself, np.transpose(np.dot(Zself, self.comaInv)))
		Aother = np.dot(Zother, np.transpose(np.dot(Zother, other.comaInv)))
		return Aself+Aother, 2*(Bself + Bother), Cself + Cother # A factor 2 was missing inB

	def getCorrectedAmplitudes(self, params):
		"""
		Get the amplitudes, shifted by zero modes according to the parameters
		"""
		if not len(params) == 2*self.nZero:
			raise IndexError("Parameter size does not match " + str(len(params)) + ' ' + str(2*self.nZero))
		amp = np.zeros((2*self.totalBins))
		for i in range(self.totalBins):
			amp[2*i  ] = self.reals[i]
			amp[2*i+1] = self.imags[i]
			for z in range(self.nZero):
				amp[2*i  ] += self.zeroModes[z][i] * params[2*z  ]
				amp[2*i+1] += self.zeroModes[z][i] * params[2*z+1]
		return amp

	def theoChi2(self, sector, funcs, params):
		"""
		Get the chi2 for the own theory (amplitude level)
		"""
		if not len(params) == 2*self.nZero + 2*len(funcs):
			raise IndexError("Parameter size does not match")
		binStart  = self.borders[sector]
		binStop   = self.borders[sector+1]
		deltas    = np.zeros((2*self.totalBins))
		count     = 0
		corrected = self.getCorrectedAmplitudes(params[:2*self.nZero])
		for bin in range(binStart, binStop):
			deltas[2*bin  ] = corrected[2*bin  ]
			deltas[2*bin+1] = corrected[2*bin+1]
			for p in range(len(funcs)):
				ampl = (params[2*(self.nZero + p)] + 1.j*params[2*(self.nZero + p)+1]) * funcs[p][bin] * self.norms[bin]**.5
				deltas[2*bin  ] -= ampl.real
				deltas[2*bin+1] -= ampl.imag
			count += 1
		return np.dot(deltas, np.dot(self.comaInv, deltas))

	def smoothnessChi2(self, other, params, useSelf = True, useOther = True):
		"""
		Gets the smoothness chi2 w.r.t. another mass bin
		"""
		if not self.sectors == other.sectors:
			raise ValueError("Sectors do not match")
		Dself  = np.zeros((2*self.totalBins))
		Dother = np.zeros((2*other.totalBins))
		paramsSelf  = params[:2*self.nZero]
		paramsOther = params[2*self.nZero:]
#		print self.nZero, other.nZero, params,"LLLLLLLLLOOPPdssa"

		ampsSelf  = self.getCorrectedAmplitudes(paramsSelf)
		ampsOther = other.getCorrectedAmplitudes(paramsOther)
		for s in range(self.nSect):
			startSelf  = self.borders[s]
			startOther = other.borders[s]
			nBins      = min(self.borders[s+1] - self.borders[s], other.borders[s+1] - other.borders[s])
			for i in range(nBins):
				Dself[2*(startSelf + i)  ] = ampsSelf[2*(i+startSelf)  ] - ampsOther[2*(i+startOther)  ]
				Dself[2*(startSelf + i)+1] = ampsSelf[2*(i+startSelf)+1] - ampsOther[2*(i+startOther)+1]

				Dother[2*(startOther + i)  ] = ampsSelf[2*(i+startSelf)  ] - ampsOther[2*(i+startOther)  ]
				Dother[2*(startOther + i)+1] = ampsSelf[2*(i+startSelf)+1] - ampsOther[2*(i+startOther)+1]
		Cself  = np.dot(Dself, np.dot(self.comaInv, Dself))
		Cother = np.dot(Dother, np.dot(other.comaInv, Dother))
		retVal = 0.
		if useSelf:
			retVal += Cself
		if useOther:
			retVal += Cother
		return retVal

	def fillTotal(self, param, hists, binRange = {}):
		if not len(hists) == self.nSect:
			raise IndexError("Number of histograms does not match the number of sectors")
		totals = self.getSectorTotals(param, binRange = binRange)
		for s in range(self.nSect):
			hists[s].SetBinContent(self.bin3pi+1, totals[s][0])
			hists[s].SetBinError(self.bin3pi+1, totals[s][1])

	def fillHistograms(self, params, hists, mode = INTENS):
		"""
		Fill coorespi=onding m3pi bin of 2D histograms
		"""
		if mode.IS_THEO and not self.hasTheo:
			print "No theory loaded, cannot fill histogram"
		if not len(hists) == self.nSect:
			raise IndexError("Histogram number mismatch")
		corrAmp = self.getCorrectedAmplitudes(params)
		for s in range(self.nSect):
			count = 0
			start = self.borders[s  ]
			stop  = self.borders[s+1]
			for i in range(start, stop):
				ampl = corrAmp[2*i] + 1.j * corrAmp[2*i+1]
				norm = self.norms[i]
				coma = np.zeros((2,2))
				jac  = np.zeros((2))
				coma[0,0] = self.coma[2*i  ,2*i  ]
				coma[0,1] = self.coma[2*i  ,2*i+1]
				coma[1,0] = self.coma[2*i+1,2*i  ]
				coma[1,1] = self.coma[2*i+1,2*i+1]
				if mode == INTENS:
					val = abs(ampl)**2
					jac[0] = 2*ampl.real
					jac[1] = 2*ampl.imag
				elif mode == INTENSNORM:
					val = abs(ampl)**2/norm
					jac[0] = 2*ampl.real/norm
					jac[1] = 2*ampl.imag/norm
				elif mode == REAL:
					val = ampl.real
					jac[0] = 1.
				elif mode == IMAG:
					val = ampl.imag
					jac[1] = 1.
				elif mode == REIMCORRELATION:
					val = coma[0,1]
				elif mode == PHASE:
					val = phase(ampl)
					if ampl.real == 0.:
						if ampl.imag > 0.:
							jac[0] = -1./ampl.imag
						else:
							jac[0] =  1./ampl.imag
					else:
						common = 1. + ampl.imag**2/ampl.real**2
						jac[0] = -ampl.imag/ampl.real**2/common
						jac[1] = 1./ampl.real/common
				elif mode == INTENSTHEO:
					val = abs(self.theo[i])**2
				elif mode == REALTHEO:
					val = self.theo[i].real
				elif mode == IMAGTHEO:
					val = self.theo[i].imag
				elif mode == PHASETHEO:
					val = phase(self.theo[i])
				else:
					raise ValueError("Unknown mode '" + mode + "'")
				err = np.dot(jac, np.dot(coma,jac))**.5
				hists[s].SetBinContent(self.bin3pi+1, count + 1, val)
				hists[s].SetBinError(self.bin3pi+1, count + 1, err)
				count += 1

	def rotateToPhaseOfBin(self, nBin):
		"""
		Rotates the global phase, that the phase of nBin is zero
		"""
		ampl = self.reals[nBin] + 1.j*self.imags[nBin]
		self.removePhase(ampl)

	def removePhase(self, amp):
		ampl = amp.real -1.j*amp.imag
		ampl/= abs(ampl)
		re = ampl.real
		im = ampl.imag
		for i in range(len(self.reals)):
			newAmpl = ampl*(self.reals[i] + 1.j*self.imags[i])
			self.reals[i] = newAmpl.real
			self.imags[i] = newAmpl.imag
#	#	CMwrite("removePhase")
		jac = np.zeros((len(self.coma), len(self.coma)))
		for i in range(len(self.coma)/2):
			jac[2*i  ,2*i  ] =   re
			jac[2*i  ,2*i+1] = - im
			jac[2*i+1,2*i  ] =   im
			jac[2*i+1,2*i+1] =   re
		self.coma         = np.dot(jac, np.dot(self.coma,np.transpose(jac)))
		self.makeComaInv()
		self.specialCOMAs = {}

	def comaIsSymmetric(self):
		"""
		Checks, if the covriance matrix is symmetric
		"""
		for i in range(2*self.totalBins):
			for j in range(2*self.totalBins):
				if not self.coma[i,j] == self.coma[j,i]:
					print i,j,self.coma[i,j],self.coma[j,i]
					return False
		return True

	def getSectorInt(self, sect):
		if isinstance(sect, int):
			if sect < 0 or sect >= len(self.sectors):
				raise ValueError("Sector index out of range: "+str(sect))
			return sect
		else:
			for s, listSect in enumerate(self.sectors):
				if listSect == sect:
					return s
		raise ValueError("Could not find sector number for '" + sect + "'")

	def toIntegerSectorMap(self, mapp):
		"""
		Changes a map of {sector : <anything>}, where the key can be string or int to a corresponding map, where the keys are ints
		"""
		if mapp == {}:
			return {}
		isInt = False
		isStr = False
		for key in mapp:
			if isinstance(key, int):
				isInt = True
			if isinstance(key,str):
				isStr = True
		if isInt and isStr:
			raise ValueError("Cannot handle both string and integer keys")
		elif isInt:
			return mapp
		elif isStr:
			newMapp = {}
			for s, sector in enumerate(self.sectors):
				if sector in mapp:
					newMapp[s] = mapp[sector]
			return newMapp
		else:
			raise ValueError("Something wrong with the mapp")


