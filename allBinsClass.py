from massBinClass import massBin
import numpy as np
import numpy.linalg as la
from modes import INTENS, PHASE


class allBins:
	def __init__(self, binStart, binStop, realHists, imagHists, normHists, indexHists, comaHistList):
		"""
		Initializes the allBins class
		"""
		self.binStart = binStart
		self.binStop  = binStop
		if not binStop - binStart == len(comaHistList):
			raise IndexError("List of COMA hists does not match")
		self.massBins = []
		count = 0
		for i in range(binStart, binStop):
			self.massBins.append(massBin(i, realHists, imagHists, normHists, indexHists, comaHistList[count]))
			count += 1
		self.chi2init = False

	def initChi2(self, sectorFuncMap):
		"""
		Initializes the chi2 function by setting the corresponting functions in all mass bins
		"""
		for mb in self.massBins:
			mb.initChi2(sectorFuncMap)
		self.chi2init = True

	def __getitem__(self, index):
		"""
		Enables looping through mass bins
		"""
		if index < 0 or index >= len(self.massBins):
			raise IndexError("Invalud index in allBins.__getitem__(" + str(index) + ")")
		return self.massBins[index]

	def __len__(self):
		"""
		Enables usage of the 'len(...)' method
		"""
		return len(self.massBins)

	def modeChi2(self, pars, mode = PHASE, mBinResolved = False):
		"""
		Evaluates a chi2 for a specific mode, other than the amplitude (PHASE or INTENS), summed over all mass bins
		"""
		if not len(pars) == self.nParAll():
			print len(pars) , self.nParAll()
			raise IndexError("Number of parameters does not match")
		if mBinResolved: 
			chi2 = []
		else:
			chi2  = 0.
		count = 0
		nPar  = self.massBins[0].nPar
		nNon  = 2*(self.massBins[0].nZero + self.massBins[0].nFunc)
		for i, mb in enumerate(self.massBins):
			par    = pars[count:count+nNon]
			if nPar > 0:
				par += pars[-nPar:]
			binC2 = mb.modeChi2(par, mode = mode)
			if mBinResolved:
				chi2.append(binC2)
			else:
				chi2 += binC2
			count += nNon
		return chi2

	def getNDF(self):
		"""
		Return NDF as list of (NDF, 2*nZero, 2*nFunc, nPar) (Simple mass bin loop)
		"""
		retVal = []
		for mb in self.massBins:
			retVal.append(mb.getNDF())
		return retVal

	def nPar(self):
		"""
		Gets the number of shape parameters (assumes to be the same in all mass bins)
		"""
		return self.massBins[0].nPar

	def nParAll(self):
		"""
		Gets the total number of parameters (Makes no asusmptions)
		"""
		nPar = 0
		for mb in self.massBins:
			nPar += 2*mb.nZero + 2*mb.nFunc
		nPar += mb.nPar
		return nPar

	def chi2(self, shapeParams = [], returnParameters = False):
		"""
		Evalueates the chi2 for the AMPL method. Non-shape parameters are calculated. If no shape parameters are given, the internally sotred ones are used.
		"""
		if not self.chi2init:
			raise RuntimeError("Chi2 not initialized, cannot evaluate")
		chi2 = 0.
		if returnParameters:
			params = []
		for mb in self.massBins:
			if returnParameters:
				c2, pars = mb.chi2(shapeParams, True)
				chi2 += c2
				params.append(pars)
			else:
				chi2 += mb.chi2(shapeParams)
		if returnParameters:
			return chi2, params
		else:
			return chi2

	def setParametersAndErrors(self, parameters, errors):
		"""
		Sets the internal parameters (and their uncertainties) 
		"""
		if not self.chi2init:
			raise RuntimeError("Chi2 not initialized, cannot set errors")
		for mb in self.massBins:
			if not mb.setParametersAndErrors(parameters, errors):
				print "Setting of error in one bin failed"
				return False
		return True

	def addZeroMode(self, borders, zeroMode, eigenvalueHist = None):
		"""
		Adds a zero-mode (Simple mass bin loop)
		"""
		for mb in self.massBins:
			mb.addZeroMode(borders, zeroMode, eigenvalueHist = None)

	def renormZeroModes(self):
		"""
		Renorms zero modes (Simple mass bin loop)
		"""
		for mb in self.massBins:
			mb.renormZeroModes()

	def removeGlobalPhaseFromComa(self):
		"""
		Removes the dorection of a gobal phase rotation from the covariance matrix (Simple mass bin loop)
		"""
		for mb in self.massBins:
			mb.removeGlobalPhaseFromComa()

	def setMassRanges(self, massRanges):
		"""
		Sets the mass ranges for the evaluation of chi2 estimators (Simple mass bin loop)
		"""
		for mb in self.massBins:
			mb.setMassRanges(massRanges)

	def nZero(self):
		"""
		Gives the total number of zero modes, summed over all mass bins
		"""
		nZero = 0
                for mb in self.massBins:
                        nZero += mb.nZero
		return nZero

	def unifyComa(self):
		"""
		Sets the covariance matrices to unity (Simple mass bin loop)
		"""
		for mb in self.massBins:
			mb.unifyComa()

	def setTheory(self, sectParsFunctMap):
		"""
		Sets internally evaluated theory curves (Simple mass bin loop). If shape paramerers are given, these are used, otherwise the internally stored ones.
		"""
		for i, mb in enumerate(self.massBins):
			spfm = {}
			for s in sectParsFunctMap:
				spfm[s] = []
				for pf in sectParsFunctMap[s]:
					p = pf[0][i]
					f = pf[1]
					if len(pf) == 2:
						spfm[s].append([p,f])
					elif len(pf) == 3:
						shapePars = pf[2][i]
						spfm[s].append([p,f,shapePars])
			mb.setTheory(spfm)

	def setTheoryFromOwnFunctions(self, parameterLists, skipZeroPars = True, restrictToRange = True):
		"""
		Sets internally evaluated theory curves from set functions
		"""
		if not self.chi2init:
			raise RuntimeError("Chi2 not inited, cannot set theory")
		if not len(parameterLists) == len(self.massBins):
			raise ValueError("Number of parameterlists does not match " + str(len(parameterLists)) + ' ' + str(len(self.massBins)))
		for i,mb in enumerate(self.massBins):
			if skipZeroPars:
				nZ = mb.nZero
				mb.setTheoryFromOwnFunctions(parameterLists[i][2*nZ:], restrictToRange = restrictToRange)
			else:
				mb.setTheoryFromOwnFunctions(parameterLists[i], restrictToRange = restrictToRange)

	def linearizeZeroModeParameters(self, pars):
		"""
		Prduces a 1D array from a 2D array of zero mode parameters
		"""
		if not len(pars) == len(self.massBins):
			raise ValueError("Number of parameters do not match")
		nZero   = self.nZero()
		linPars = np.zeros((2*nZero))
		count   = 0
		for i,mb in enumerate(self.massBins):
			nZ = mb.nZero
			for z in range(nZ):
				linPars[2*(z+count)  ] = pars[i][2*z  ]
 				linPars[2*(z+count)+1] = pars[i][2*z+1]
			count += nZ
		return linPars

	def getSmoothnessABC(self):
		"""
		Retuns the A, B and C of the smoothness method, wehere: chi2(p) = pAp + Bp + C
		"""
		nZero = self.nZero()
		C = 0.
		B = np.zeros((2*nZero))
		A = np.zeros((2*nZero, 2*nZero))
		countZero = 0
		for i in range(self.binStop - self.binStart - 1 ):
			a,b,c = self.massBins[i].getSmoothnessABC(self.massBins[i+1])
			C += c
			for j in range(self.massBins[i+1].nZero + self.massBins[i].nZero):
				B[2*(countZero + j)  ] += b[2*j  ]
				B[2*(countZero + j)+1] += b[2*j+1]
				for k in range(self.massBins[i+1].nZero + self.massBins[i].nZero):
					A[2*(countZero + j)  ,2*(countZero + k)  ] += a[2*j  ,2*k  ]
                                        A[2*(countZero + j)  ,2*(countZero + k)+1] += a[2*j  ,2*k+1]
                                        A[2*(countZero + j)+1,2*(countZero + k)  ] += a[2*j+1,2*k  ]
                                        A[2*(countZero + j)+1,2*(countZero + k)+1] += a[2*j+1,2*k+1]
			countZero += self.massBins[i].nZero
#		val,vec = la.eig(A)
#		for v in val:
#			print v,
#			if v.real  < 0.:
#				print "!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#			else:
#				print
		return A,B,C

	def writeZeroModeCoefficients(self, coefficients, outFileName, tBin = "<tBin>"):
		"""
		Writes the zero mode coefficients to a text file (Simple mass bin loop)
		"""
		nZero = self.nZero()
		cmplx = False
		if len(coefficients) == nZero:
			cmplx = True
		elif not len(coefficients) == 2*nZero:
			raise IndexError("Number of coefficients for "+str(self.nZero)+" modes does not match (neither real nor complex): " +str(len(coefficients)))
		tBin  = str(tBin)		
		count = 0
		for mb in self.massBins:
			nn = mb.nZero
			if not cmplx:
				nn *= 2
			mb.writeZeroModeCoefficients(coefficients[count:count+nn], outFileName, tBin, 'a')
			count += nn

	def printNzero(self):
		"""
		Prins the numbers of zero modes in all bins
		"""
		nZeros = []
		for mb in self.massBins:
			nZeros.append(mb.nZero)
		print nZeros

	def removeZeroModeFromComa(self):
		"""
		Removes the direction of a shift in the zero mode direction from the covariance matrix (Simple mass bin loop)
		"""
		for mb in self.massBins:
			mb.removeZeroModeFromComa()

	def smothnessChi2(self,params):
		"""
		Gives the chi2 for the smoothness method
		"""
		countZero = 0
		chi2 = 0.
		for i in range(self.binStop - self.binStart - 1 ):
			chi2 += self.massBins[i].smothnessChi2(self.massBins[i+1], params[2*countZero:2*(countZero + self.massBins[i+1].nZero + self.massBins[i].nZero)])
			countZero += self.massBins[i].nZero
		return chi2

	def resolvedSmoothnessChi2(self, params):
		"""
		Gives the chi2 for the smoothness method, resolved by m bins (Contribution of every bin's COMA)
		"""	
		countZero = 0
		chi2s     = [0.] * (self.binStop - self.binStart)
		for i in range(self.binStop - self.binStart):
			if not i == 0:
				chi2s[i]  += self.massBins[i-1].smoothnessChi2(self.massBins[i], params[2*countZero:2*(countZero + self.massBins[i].nZero + self.massBins[i-1].nZero)], useSelf = False)
				countZero += self.massBins[i].nZero
			if not i == self.binStop - self.binStart-1:
				chi2s[i] += self.massBins[i].smoothnessChi2(self.massBins[i+1], params[2*countZero:2*(countZero + self.massBins[i+1].nZero + self.massBins[i].nZero)], useOther = False)
		return chi2s


	def rotateToPhaseOfBin(self, nBin):
		"""
		Rotates the global phase, that the phase of nBin is zero (Simple mass bin loop)
		"""
		for mb in self.massBins:
			mb.rotateToPhaseOfBin(nBin)

	def fillHistograms(self, params, hists, mode = INTENS):
		"""
		Fills 2D histograms
		"""
		parCount = 0
		for mb in self.massBins:
			nz = 2*mb.nZero
			mb.fillHistograms(params[parCount:parCount+nz], hists, mode = mode)
			parCount += nz		
	
