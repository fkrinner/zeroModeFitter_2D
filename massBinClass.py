import pyRootPwa
import numpy as np
import numpy.linalg as la
import utils
from utils import getZeroHistSectors, normVector, getZeroModeNumber
from cmath import phase, pi
from modes import INTENS, PHASE, REAL, IMAG, INTENSNORM, INTENSTHEO, REALTHEO, IMAGTHEO, PHASETHEO, REIMCORRELATION

class massBin:
	def __init__(self, nBin, realHists, imagHists, normHists, indexHists, coma):
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
		self.coma       = np.zeros((2*self.totalBins,2*self.totalBins))
		self.binCenters = np.zeros((self.totalBins))
		self.numLim     = 1.e-11
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
						count2 += 1
				count +=1
		self.makeComaInv()
		self.borders = [0]
		for i in range(self.nSect):
			self.borders.append(self.borders[-1] + self.nBins[i])
		self.nZero           =  0
		self.zeroModes       = [ ]
		self.zeroModeNumbers = [ ]
		self.zeroModeTitles  = [ ]
		self.zeroEigenvalues = [ ]
		self.hasTheo         = False
		self.hasMassRange    = False
		self.chi2init        = False
		self.specialCOMAs    = { }

	def initChi2(self, sectorFuncMap):
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
		self.nFunc = 0
		for val in self.nFuncs:
			self.nFunc += val
		self.chi2init = True

	def chi2(self, pars = [], returnParameters = False):
		if not self.chi2init:
			raise RuntimeError("chi2 not inited, cannot evaluate")
		if len(pars) > 0:
			self.setShapeParameters(pars)
		A,B,C = self.getOwnTheoryABC()
		pars  = -np.dot(B, la.pinv(A + np.transpose(A)))
		chi2  = np.dot(pars,np.dot(A,pars)) + np.dot(pars,B) + C
		if returnParameters:
			return chi2, pars
		else:
			return chi2

	def nParAll(self):
		return 2*self.nZero + 2*self.nFunc + self.nPar

	def phaseChi2(self, pars):
		return self.modeChi2(pars, PHASE)

	def intensChi2(self, pars):
		return self.modeChi2(pars, INTENS)

	def modeChi2(self, pars, mode = PHASE):
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
				elif mode == INTENS:
					deltas[bin] = corrAmpl[2*bin]**2 + corrAmpl[2*bin+1]**2 - abs(self.theo[bin])**2


		coma = utils.pinv(self.getSpecialComa(mode), self.numLim)
		retVal = np.dot(deltas, np.dot(coma, deltas))

		removePhaseDirection = True
		if removePhaseDirection:
			entry          = (1./len(deltas))**.5 # The direction of a global phase rotation is (1., 1., 1., 1., ...) / len(deltas)**.5
			sp             = 0.
			for i in range(len(deltas)):
				sp += entry * deltas[i]
			retVal += sp**2*retVal

		return retVal

	def getSpecialComa(self, mode = PHASE):
		if mode in self.specialCOMAs:
			return self.specialCOMAs[mode]
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
		retVal = np.dot(jacobian, np.dot(self.coma, np.transpose(jacobian)))
		self.specialCOMAs[mode] = retVal
		return retVal

	def setShapeParameters(self, pars):
		if not self.chi2init:
			raise RuntimeError("chi2 not inited, no knowledge about shape parameters")
		if not len(pars) == self.nPar:
			raise ValueError("Number of shape parameters does not match len(pars) = " + str(len(pars)) + " != self.nPar = " + str(self.nPar))
		countPar = 0
		for s in range(self.nSect):
			if not s in self.funcs:
				continue
			for f in self.funcs[s]:
				f.setParameters(pars[countPar:countPar + f.nPar])
				countPar += f.nPar

	def setParametersAndErrors(self, pars, errors):
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
		self.theo = np.zeros((self.totalBins), dtype = complex)

	def setTheory(self, sectParFunctMap = {}):
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
		for i in range(len(self.coma)):
			for j in range(len(self.coma)):
				if i == j:
					self.coma[i,j] = 1.
					self.comaInv[i,j] = 1.
				else:
					self.coma[i,j] = 0.
					self.comaInv[i,j] = 0.
		self.specialCOMAs = {}

	def makeComaInv(self):
		self.comaInv = utils.pinv(self.coma, self.numLim)


	def removeZeroModeFromComa(self):
		if len(self.zeroModes) == 0:
			return
		dim = len(self.coma)
		transformationMatrix = np.identity(dim)
		for z in range(self.nZero):
			for i in range(dim/2):
				for j in range(dim/2):
					transformationMatrix[2*i  ,2*j  ] -= self.zeroModes[z][i] * self.zeroModes[z][j]
					transformationMatrix[2*i+1,2*j+1] -= self.zeroModes[z][i] * self.zeroModes[z][j]
		self.coma         = np.dot(transformationMatrix, np.dot(self.coma, transformationMatrix)) # no transpose needed, since transformationMatrix is symmetric
		self.makeComaInv()
		self.specialCOMAs = {}

	def getGlobalPhaseDirection(self):
		retVal = np.zeros((2*self.totalBins))
		for i in range(self.totalBins):
			retVal[2*i  ] = -self.imags[i]
			retVal[2*i+1] =  self.reals[i]

		return normVector(retVal)

	def removeGlobalPhaseFromComa(self):
		phaseDirection = self.getGlobalPhaseDirection()
		dim = 2*self.totalBins
		transformationMatrix = np.identity(dim)
		for i in range(dim):
			for j in range(dim):
				transformationMatrix[i,j] -= phaseDirection[i] * phaseDirection[j]
		self.coma         = np.dot(transformationMatrix, np.dot(self.coma, transformationMatrix)) # no transpose needed, since transformationMatrix is symmetric
		self.makeComaInv()
		self.specialCOMAs = {}

	def setMassRanges(self, massRanges):
		self.massRanges   = self.toIntegerSectorMap(massRanges)
		self.hasMassRange = True

	def setTheoryFromOwnFunctions(self, parameterList, restrictToRange = True):
		if not self.chi2init:
			raise RuntimeError("Chi2 not inited, does not have own functions")
		if not len(parameterList) == 2*self.nFunc:
			raise ValueError("Number of parameters does not match")
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

	def writeAmplitudeAndArgandFile(self, params, outFileNameBase):
		pass

		

	def writeZeroModeCoefficients(self, coefficients, outFileName, tBin = "<tBin>", mode = 'a'):
		tBin  = str(tBin)
		cmplx = False
		if len(coefficients) == self.nZero:
			cmplx = True
		elif not len(coefficients) == 2*self.nZero:
			raise IndexError("Number of coefficients for "+str(self.nZero)+" modes does not match (netiher real nor complex): " +str(len(coefficients)))
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
		for z in range(self.nZero):
			for i in range(self.totalBins):
				self.zeroModes[z][i] *= self.norms[i]**.5

	def addZeroMode(self, modeBorders, modeHist, eigenvalueHist = None):
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

	def getTheoryABC(self, sector, parametrizations, massRange = None):
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
		if not len(params) == 2*self.nZero:
			raise IndexError("Parameter size does not match")
		amp = np.zeros((2*self.totalBins))
		for i in range(self.totalBins):
			amp[2*i  ] = self.reals[i]
			amp[2*i+1] = self.imags[i]
			for z in range(self.nZero):
				amp[2*i  ] += self.zeroModes[z][i] * params[2*z  ]
				amp[2*i+1] += self.zeroModes[z][i] * params[2*z+1]
		return amp

	def theoChi2(self, sector, funcs, params):
		if not len(params) == 2*self.nZero + 2*len(funcs):
			raise IndexError("Parameter size does not match")
		binStart = self.borders[sector]
		binStop  = self.borders[sector+1]
		deltas   = np.zeros((2*self.totalBins))
		count = 0
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

	def smothnessChi2(self, other, params):
		if not self.sectors == other.sectors:
			raise ValueError("Sectors do not match")
		Dself  = np.zeros((2*self.totalBins))
		Dother = np.zeros((2*other.totalBins))
		paramsSelf  = params[:2*self.nZero]
		paramsOther = params[2*self.nZero:]
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
		return Cself + Cother

	def fillHistograms(self, params, hists, mode = INTENS):
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
		ampl = self.reals[nBin] - 1.j*self.imags[nBin]
		ampl/= abs(ampl)
		re = ampl.real
		im = ampl.imag
		for i in range(len(self.reals)):
			newAmpl = ampl*(self.reals[i] + 1.j*self.imags[i])
			self.reals[i] = newAmpl.real
			self.imags[i] = newAmpl.imag

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
		for i in range(2*self.totalBins):
			for j in range(2*self.totalBins):
				if not self.coma[i,j] == self.coma[j,i]:
					print i,j,self.coma[i,j],self.coma[j,i]
					return False
		return True

	def toIntegerSectorMap(self, mapp):
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


