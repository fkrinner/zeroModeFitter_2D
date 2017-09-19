from rootfabi import root_open, GetKeyNames
import numpy as np
import numpy.linalg as la
from allBinsClass import allBins
from modes import PHASE, AMPL, SMOOTH, NONE, REAL, IMAG, REIMCORRELATION, INTENSTHEO, REALTHEO, IMAGTHEO, PHASETHEO
import os, sys
from utils import zeroForSectors, getZeroHistBorders, sumUp, INF, renormToBinWidth, divideAbyB, get3PiHistogram
import parameterTrackingParameterizations as ptc
import scipy.optimize
from removeZeroModes import removeCertainZeroModes
from tBinHolder import tBinHolder
from random import random
from LaTeX_strings import getProperWaveName, getProperDataSet
from resultViewerClass import resultViewer
from utils import loadAmplsTM, changeReferenceWave
from estimateErrors import estimateErrors,estimateErrors2,estimateErrors3


			

class amplitudeAnalysis:
	"""
	Class holding all information to fix zero mode coefficients in several ways
	and do this consistently
	"""
	def __init__(self, inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = {}, zeroFileName = ""):
		"""
		Simple contructor
		"""
		self.inFileName = inFileName
		if zeroFileName == "":
			self.zeroFileName = inFileName
		else:
			self.zeroFileName = zeroFileName
		if not os.path.isfile(inFileName):
			raise IOError("File '" + inFileName + "' does not exist")
		self.sectors        = sectors
		self.sectorRangeMap = sectorRangeMap
		self.waveModel      = waveModel
		self.startBin       = startBin
		self.stopBin        = stopBin
		self.tBins          = tBins
		self.model          = tBinHolder()
		self.flags          = {}
		self.mode           = NONE

	def IS(self, flag):
		"""
		Checks is a certain flag is set and still True
		"""
		if not flag in self.flags:
			return False
		return self.flags[flag]

	def SET(self, flag):
		"""
		Sets a certain flag
		"""
		self.flags[flag] = True

	def fitShapeParametersForBinRange(self, startPars, tBinsToEvaluate, mBinsToEvaluate, zeroModeParameters = []):
		if len(zeroModeParameters) > 0:
			self.setZeroModeParameters(zeroModeParameters)
		self.model.setBinsToEvalueate(tBinsToEvaluate,mBinsToEvaluate)
		pars = startPars[:]
		res       = scipy.optimize.minimize(self.model.fixedZMPchi2, pars)
		hi        = res.hess_inv
		errs      = []
		startErrs = []
		pars = res.x[:]
		for i in range(len(res.x)):
#			print hi[i,i],"::::::::::::::::::::{{}"
			errs.append((2.*hi[i,i])**.5)
		errs = estimateErrors2(self.model.chi2, pars, errs)
		return res.x, errs


	def setZeroModeParameters(self,params):
		self.model.setZeroModeParameters(params)

	def loadData(self, loadIntegrals = False, phaseFile = "", referenceWave = ""):
		"""
		Loads the rquired data from a ROOT file
		"""
		if not phaseFile == "":
			raise RuntimeError("The method with the phase file is outdated... use the referenceWave argument")
#			phaseAmplitudes = loadAmplsTM(phaseFile)

		for tBin in self.tBins:
			with root_open(self.inFileName, "READ") as inFile:
				histNames = GetKeyNames(inFile)
				histListReal = []
				histListImag = []
				histListNorm = []
				histListIndx = []
				for sector in self.sectors:
					realHistName = sector+"_"+str(tBin)+"_real"
					histReal = inFile.Get(realHistName)
					if not histReal:
						raise IOError("Could not get '" + realHistName + "' from '" + self.inFileName + "'")
					histListReal.append(histReal)

					imagHistName = sector+"_"+str(tBin)+"_imag"
					histImag = inFile.Get(imagHistName)
					if not histImag:
						raise IOError("Could not get '" + imagHistName + "' from '" + self.inFileName + "'")
					histListImag.append(histImag)

					normHistName = sector+"_"+str(tBin)+"_norm"
					histNorm = inFile.Get(normHistName)
					if not histNorm:
						raise IOError("Could not get '" + normHistName + "' from '" + self.inFileName + "'")
					histListNorm.append(histNorm)

					indexHistName = sector+"_"+str(tBin)+"_index"
					histIndx = inFile.Get(indexHistName)
					if not histIndx:
						raise IOError("Could not get '" + indexHistName + "' from '" + self.inFileName + "'")
					histListIndx.append(histIndx)
				histsToFill  = []
				comaHists    = []
				intHistsReal = []
				intHistsImag = []
				for mBin in range(self.startBin, self.stopBin):
					comaHistName = "COMA_"+str(tBin)+"_"+str(mBin)
					comaHist = inFile.Get(comaHistName)
					if not comaHist:
						raise IOError("Could not get '" + comaHistName + "' from '" + self.inFileName + "'")
					comaHists.append(comaHist)
					if loadIntegrals:
						realHistName = "INTEGRAL_r_"+str(tBin)+"_"+str(mBin)
						realHist = inFile.Get(realHistName)
						if not realHist:
							raise IOError("Could not get '" + realHistName + "' from '" + self.inFileName + "'")
						intHistsReal.append(realHist)
						imagHistName = "INTEGRAL_i_"+str(tBin)+"_"+str(mBin)
						imagHist = inFile.Get(imagHistName)
						if not imagHist:
							raise IOError("Could not get '" + imagHistName + "' from '" + self.inFileName + "'")
						intHistsImag.append(imagHist)
				if not referenceWave == "":
					refHistReal  = inFile.Get(referenceWave+"_"+str(tBin)+"_real")
					if not refHistReal:
						raise IOError("Could not get '" + referenceWave+"_"+str(tBin)+"_real" + "' from '" + self.inFileName + "'")
					refHistImag  = inFile.Get(referenceWave+"_"+str(tBin)+"_imag")
					if not refHistImag:
						raise IOError("Could not get '" + referenceWave+"_"+str(tBin)+"_imag" + "' from '" + self.inFileName + "'")
					refHistIndex = inFile.Get(referenceWave+"_"+str(tBin)+"_index")
					if not refHistIndex:
						raise IOError("Could not get '" + referenceWave+"_"+str(tBin)+"_index" + "' from '" + self.inFileName + "'")
					changeReferenceWave(histListReal, histListImag, histListIndx, comaHists, refHistReal, refHistImag, refHistIndex, self.startBin, self.stopBin)


				ab = allBins(self.startBin, self.stopBin, histListReal, histListImag, histListNorm, histListIndx, comaHists,intHistsReal,intHistsImag)
				for h in histListReal:
					h.SetDirectory(0)
				self.histListReal = histListReal

			with root_open(self.zeroFileName, "READ") as inFile:
				zeroCount = 0
				zeroHistList  = []
				eigenHistList = []
				while True:
					zeroName  = "zero"+str(zeroCount)+"_"+str(tBin)
					eigenName = "eigen"+str(zeroCount)+"_"+str(tBin)
					zeroHist = inFile.Get(zeroName)
					if not zeroHist:
						break
#					print "Adding zero-mode"
					zeroCount += 1
					if not zeroForSectors(self.sectors, zeroHist.GetTitle()):
						continue
					zeroHistList.append(zeroHist)

					eigenHist = inFile.Get(eigenName)
					if eigenHist:
						eigenHistList.append(eigenHist)
				if (not len(eigenHistList) == 0) and (not len(eigenHistList) == len(zeroHistList)):
					raise ValueError("Number of eigenvalue histograms does not match, but is also nonzero")
				removeCertainZeroModes(zeroHistList, eigenHistList)
				for zeroHist in zeroHistList:
					borders = getZeroHistBorders(zeroHist)
					ab.addZeroMode(borders, zeroHist)
				if not phaseFile == "":
					ab.removePhases(phaseAmplitudes[tBin][self.startBin:self.stopBin])
				self.model.addBin(ab)
			self.SET('loaded')

	def writeZeroModeCoefficients(self, fileName, coeffs, tIn):
		if not len(self.model) == len(coeffs):
			raise ValueError("Dimension mismatch on t level")
		for t,tBin in enumerate(self.model):
			ccc = []
			for p in coeffs[t]:
				for v in p:
					ccc.append(v)
			tBin.writeZeroModeCoefficients(ccc, fileName, str(tIn))

	def getTotalHists(self, params, binRange = {}):
		hists = []
		for t, tBin in enumerate(self.model):
			tLine = []
			for m in range(len(self.sectors)):
				tLine.append(get3PiHistogram(self.sectors[m] + "_t"+ str(self.tBins[t])))
			hists.append(tLine)
		self.fillTotal(params, hists, binRange = binRange)
		return hists


	def fillTotal(self, params, hists, binRange = {}):
		for t,tBin in enumerate(self.model):
			linPar = []
			for m in params[t]:
				for v in m:
					linPar.append(v)
			tBin.fillTotal(linPar, hists[t], binRange = binRange)

	def getZeroModeSignature(self):
		"""
		Gets the number of zero modes in every t and m bin 
		"""
		retVal = []
		for tBin in self.model:
			tVals = []
			for mBin in tBin:
				tVals.append(mBin.nZero)
			retVal.append(tVals)
		return retVal

	def setZeroModeSignature(self, signature, position):
		ownSig = self.getZeroModeSignature()
		if not len(ownSig) == len(signature):
			raise ValueError("Dimension mismatch at t' level")
		for t in range(len(signature)):
			if not len(ownSig[t]) == len(signature[t]):
				raise ValueError("Dimension mismatch at m level")
			for m in range(len(signature[t])):
				if ownSig[t][m] > signature[t][m]:
					raise ValueError("Cannot reduce the number od zero-modes")
		self.SET("artificial_signature")
		self.artSig = signature
		self.additionalPosition = position

	def reduceArtificialSignatureToOwn(self, zeroModeParams):
		if not self.IS("artificial_signature"):
			return zeroModeParams
		ownSig = self.getZeroModeSignature()
		artSig = self.artSig
		a      = self.additionalPosition
		retVal = []
		if not len(ownSig) == len(zeroModeParams):
			raise ValueError("Dimension mismatch at t' level")
		for t in range(len(ownSig)):
			if not len(ownSig[t]) == len(zeroModeParams[t]):
				raise ValueError("Dimension mismatch at m level")
			tLine = []
			for m in range(len(ownSig[t])):
				if ownSig[t][m] == artSig[t][m]:
					tLine.append(zeroModeParams[t][m])
				else:
					mLine = np.zeros((2*ownSig[t][m]))
					for i in range(2*a):
						mLine[i] = zeroModeParams[t][m][i]
					for i in range(2*a, 2*ownSig[t][m]):
						mLine[i] = zeroModeParams[t][m][i+2]
					tLine.append(mLine)
			retVal.append(tLine)
		return retVal

	def createArtificialSignature(self, zeroModeParams):
		if not self.IS("artificial_signature"):
			return zeroModeParams
		ownSig = self.getZeroModeSignature()
		artSig = self.artSig
		a      = self.additionalPosition
		retVal = []
		if not len(ownSig) == len(zeroModeParams):
			raise ValueError("Dimension mismatch at t' level")
		for t in range(len(ownSig)):
			if not len(ownSig[t]) == len(zeroModeParams[t]):
				raise ValueError("Dimension mismatch at m level")
			tLine = []
			for m in range(len(ownSig[t])):
				if ownSig[t][m] == artSig[t][m]:
					tLine.append(zeroModeParams[t][m])
				else:
					mLine = np.zeros((2*artSig[t][m]))
					for i in range(2*a):
						mLine[i] = zeroModeParams[t][m][i]
					mLine[2*a]   = 0.
					mLine[2*a+1] = 0.
					for i in range(2*a, 2*ownSig[t][m]):
						mLine[i+2] = zeroModeParams[t][m][i]
					tLine.append(mLine)
			retVal.append(tLine)
		return retVal


	def finishModelSetup(self):
		"""
		Finishes model setup
		"""
		if not self.IS('loaded'):
			raise RuntimeError("Data not loaded, call 'loadData(...)' first")
		for tBin in self.model:

#			tBin.rotateToPhaseOfBin(10)
#			tBin.removeZeroModeFromComa()
#			tBin.unifyComa()

			tBin.initChi2(self.waveModel)
			tBin.setMassRanges(self.sectorRangeMap)
		self.SET('setup')

	def  compareTwoZeroModeCorrections(self, params1, params2):
		"""
		Compares two different corrections of the zero-modes
		"""
		print "--in analysis class"
		print '--',params1
		print '--',params2
		print "--out analysis class"
		return self.model.compareTwoZeroModeCorrections(params1, params2)

	def removeGlobalPhaseFromComa(self):
		"""
		Projects the direction of a global phase rotation from the covatiance matrix
		"""
		if not self.IS('setup'):
			raise RuntimeError("Model not set up, call 'finishModelSetup()' first")
		for tBin in self.model:
			tBin.removeGlobalPhaseFromComa()
		self.SET('noPhaseDir')

	def unifyComa(self):
		"""
		Sets covariance matrices to unity
		"""
		if not self.IS("setup"):
			raise RuntimeError("Model not set up, call 'finishModelSetup()' first")
		for tBin in self.model:
			tBin.unifyComa()
		self.SET("UnitCOMA")

	def fitShapeParameters(self, printResult = True):
		"""
		Performs a fit to determine the shape parameters in the AMPL mode
		"""
		if not self.IS('setup'):
			raise RuntimeError("Model has not been set up, call 'finishModelSetup()' first.")
		totalPars = []
		classParameters = []
		for k in self.waveModel:
			for f in self.waveModel[k]:
				totalPars       += f.getParameters()
				classParameters += f.parameters
		if not len(totalPars) == 0:
#			print "totalPars", totalPars			

#			delta = 0.01
#			print self.model.chi2(totalPars),":PP"
#			totalPars[0] += delta
#			print self.model.chi2(totalPars),":PP"
#			totalPars[0] -= delta
#			totalPars[1] += delta
#			print self.model.chi2(totalPars),":PP"
#			totalPars[1] -= delta
			res       = scipy.optimize.minimize(self.model.chi2, totalPars)
			hi        = res.hess_inv
			errs      = []
			startErrs = []
			pars = res.x[:]
			for i in range(len(res.x)):
				errs.append((2.*hi[i,i])**.5)
			print "errs before",errs
			print "do not estimate errors!!!"
#			errs = estimateErrors3(self.model.chi2, pars, errs)
			for i in range(len(res.x)):
				print "-/-/-/-/-/-/-/-/-/-",i,"-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-"

				pars[i] += errs[i]
				print self.model.chi2(pars) - res.fun
				pars[i] -= 2*errs[i]
				print self.model.chi2(pars) - res.fun
				pars[i] += errs[i]
			print "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-"

			
			if res.x[0] < 1.05:
				print "Compare with rho parameters"
				print "sigM =", (abs(res.x[0]-0.7690)/errs[0])
				print "sigG =", (abs(res.x[1]-0.1509)/errs[1])
			else:
				print "Compare with f2 parameters"
				print "sigM =", (abs(res.x[0]-1.2751)/errs[0])
				print "sigG =", (abs(res.x[1]-0.1851)/errs[1])


			       

#			estedErrs = estimateErrors(self.model.chi2, res.x, startErrs)
#			print "---------------------------------------"
#			print errs
#			print estedErrs
#			print "---------------------------------------"

			if not self.model[0].setParametersAndErrors(res.x, errs):
				raise RuntimeError("Could not set parameter errors")
			print "The final Chi2 is:",res.fun
			if printResult:
				for parameter in classParameters:
					print parameter
			self.chi2 = res.fun
			self.SET('hasFitResult')
			self.fitParameters = res.x[:]
		else:
			self.SET('hasFitResult')
			self.fitParameters = []
			self.chi2 = self.model.chi2(self.fitParameters)

	def calculateNonShapeParameters(self):
		"""
		Calculates the non-shape parameters---i.e. the zero-mode and model coefficients---in the AMPL mode
		"""
		if not self.IS('hasFitResult'):
			raise RuntimeError("No fit parameters, call 'fitShapeParameters(...)' first")
		self.nonShapeParameters = []
		for tBin in self.model:
			chi2, params = tBin.chi2(returnParameters = True)
			self.nonShapeParameters.append(params)
		self.SET('hasNonShape')

	def removeZeroModeFromComa(self):
		"""
		Simple loop over model bins. Removes the zero-modes from the covariance matrix
		"""
		for tt in self.model:
			tt.removeZeroModeFromComa()

	def getZeroModeParameters(self):
		"""
		Returns the zero mode parameters for the AMPL mode
		"""
		if not self.IS('hasNonShape'):
			raise RuntimeError("No non-shape parameters calculated, call 'calculateNonShapeParameters()' first")
		zeroModePars = []
		for t,tBin in enumerate(self.model):
			mBinPars = []
			for m,mBin in enumerate(tBin):
				params = self.nonShapeParameters[t][m][:2*mBin.nZero]
				oneBinParams = np.zeros((len(params)))
				for i,p in enumerate(params):
					oneBinParams[i] = p
				mBinPars.append(oneBinParams)
			zeroModePars.append(mBinPars)
		return zeroModePars

	def getSmoothnessZeroModeParameters(self):
		"""
		Returns the zero mode parameters for the smoothness method
		"""
		if not self.IS('setup'):
			raise RuntimeError("Mode has not been set up, call, 'finishModelSetup()' first.")
		zeroModePars = []
		for t, tBin in enumerate(self.model):
			a,b,c = tBin.getSmoothnessABC()
			if len(b) == 0:
				pars = np.zeros((0))
			else:
				pars  = -1. * np.dot(la.pinv(a + np.transpose(a)), b)
			bins  = []
			count = 0
			for mBin in tBin:
				nZ = mBin.nZero
				pp = np.zeros((2*nZ))
				for i in range(2*nZ):
					pp[i] = pars[count + i]
				count += 2*nZ
				bins.append(pp)
			zeroModePars.append(bins)
		return zeroModePars

	def getSmoothnessChi2s(self, pars):
		"""
		Evaluates the chi2 of the smoothness method, returns t' bin resolved
		"""
		if not len(pars) == len(self.model):
			raise IndexError("Number of t'-bins does not match")
		chi2s = []
		for t,tBin in enumerate(self.model):
			if not len(tBin) == len(pars[t]):
				raise IndexError("Number of m-bins does not match at t'-bin #" + str(t))
			a,b,c = tBin.getSmoothnessABC()
			coeff  = np.zeros((len(b)))
			count = 0
			for mm in pars[t]:
				for i in range(len(mm)):
					coeff[count] = mm[i]
					count += 1
			if not count == len(coeff):
				raise IndexError("Index mismatch found")
			chi2 = np.dot(coeff, np.dot(a,coeff)) + np.dot(coeff, b) + c
			chi2s.append(chi2)
		return chi2s

	def resolvedSmoothnessEvaluation(self, pars):
		"""
		Evaluates the chi2 of the smoothness method, returns t' amd m bin resolved
		"""
		if not len(pars) == len(self.model):
			raise IndexError("Number of t' bins does not match")
		chi2s = []
		for t, tBin in enumerate(self.model):
			linPars = []
			for mm in pars[t]:
				for v in mm:
					linPars.append(v)
			if not len(tBin) == len(pars[t]):
				raise IndexError("Number of m-bins does not match at t'-bin #" + str(t))
			chi2s.append(tBin.resolvedSmoothnessChi2(linPars))
		return chi2s

	def getPhaseChi2s(self, params, mBinResolved = False):
		"""
		Gives the chi2 of the PHASE mode for a set of zero-mode parameters
		"""
		if not self.IS('hasNonShape'):
			raise RuntimeError("Non-shape parameters not calculated, call 'phaseFit()', first'")
		wholeParams = []
		for t, tBinParameters in enumerate(self.nonShapeParameters):
			for m, mBinParameters in enumerate(tBinParameters):
				nZero = self.model[t][m].nZero
				for i,v in enumerate(mBinParameters):
					if i < 2*nZero:
						wholeParams.append(params[t][m][i])
					else:
						wholeParams.append(v)

		for v in self.fitParameters:
			wholeParams.append(v)
		return self.model.modeChi2(np.asarray(wholeParams), tBinResolved = True, mBinResolved = mBinResolved)

	def singleBinPhaseFit(self, nAttempts = 1):
		"""
		Performs a bin-wise phase fit, if no global shape parameter are present, to be quicker than 'phaseFit()' (Fits shape and non-shape parameters)
		"""
		if not self.IS('setup'):
			raise RuntimeError("Model has not been set up, call 'finishModelSetup()' first.")
		if not self.model.nPar() == 0:
			raise RuntimeError("singleBinPhaseFit() cannot be called using a model with free parameters")
		totalChi2 = 0.
		self.nonShapeParameters = []
		for t,tBin in enumerate(self.model):
			tLine = []
			for m,mBin in enumerate(tBin):
				print "tBin: "+str(t)+ "; mBin: "+str(m)
				bestChi2   = INF
				bestAttepmt = 0
				bestPars   = []
				for a in range(nAttempts):
					startPar = [2.*random() - 1. for _ in range(mBin.nParAll())]
					res      = scipy.optimize.minimize(mBin.modeChi2, startPar)
					if res.fun < bestChi2:
						bestChi2    = res.fun
						bestPar     = res.x[:]
						bestAttepmt = a
				print "adding to phase chi2:", bestChi2
				totalChi2 += bestChi2
				tLine.append(bestPar)
				print "The best attempt was attempt #" + str(bestAttepmt)
			self.nonShapeParameters.append(tLine)
		self.fitParameters = []
		print "The final Chi2 in bin by bin minimization was:",totalChi2
		self.SET("hasFitResult")
		self.SET("hasNonShape")

	def phaseFit(self, printResult = True):
		"""
		Performs the specified fit to the phase of the model (Fits shape and non-shape parameters)
		"""
		if not self.IS('setup'):
			raise RuntimeError("Model has not been set up, call 'finishModelSetup()' first.")
		if self.model.nPar() == 0:
			print "The model does not have any parameters, switch to binwise fit"
			self.singleBinPhaseFit(nAttempts = 5)
			return
		totalPars = [2.*random()-1 for _ in range(self.model.nParAll() - self.model.nPar())]
		classParameters = []
		for k in self.waveModel:
			for f in self.waveModel[k]:
#				print "========================================"
#				print type(totalPars)
#				print type(f.getParameters())
#				print "========================================"
				totalPars       += list(f.getParameters())
				classParameters += f.parameters
		res    = scipy.optimize.minimize(self.model.modeChi2, totalPars)
		hi     = res.hess_inv
		errs   = []
		for i in range(len(res.x)):
			errs.append((2*hi[i,i])**.5)
		if -self.model.nPar() > 0:
			if not self.model[0].setParametersAndErrors(res.x[-self.model.nPar():], errs[-self.model.nPar():]):
				raise RuntimeError("Could not set parameter errors")
		print "The final Chi2 is:",res.fun
		if printResult:
			for parameter in classParameters:
				print parameter
		self.chi2 = res.fun
		self.SET('hasFitResult')
		if self.model.nPar() > 0:
			self.fitParameters = res.x[-self.model.nPar():]
		else:
			self.fitParameters = []
		count = 0
		self.nonShapeParameters = []
		for tt in self.model:
			tLine = []
			for mm in tt:
				nn = (mm.nZero + mm.nFunc)*2
				tLine.append(res.x[count:count+nn])
				count += nn
			self.nonShapeParameters.append(tLine)
		self.SET('hasNonShape')

	def evaluateZeroModeParameters(self, params):
		"""
		Evaluates the chi2 for a set of zero mode parameters (this is AMPL mode)
		"""
		if not self.IS('hasFitResult'):
			raise RuntimeError("No fit parameters, call 'fitShapeParameters(...)' first")
		self.model.chi2(self.fitParameters)
		if not len(params) == len(self.model):
			print len(params) , len(self.model)
			raise IndexError("Number of tBins does not match")
		chi2s = []
		for t, tBin in enumerate(self.model):
			if not len(params[t]) == len(tBin):
				raise IndexError("Number of mBins does not match")
			c2s = []
			for m,mBin in enumerate(tBin):
				if not len(params[t][m]) == 2*mBin.nZero:
					raise IndexError("Number of zeroModes does not match " + str(len(params[t][m]))+  " vs. " + str(2*mBin.nZero))
				wholeParams = np.copy(self.nonShapeParameters[t][m])
				for i in range(2*mBin.nZero):
					wholeParams[i] = params[t][m][i]
				a,b,c = mBin.getOwnTheoryABC()
				for i in range(len(b)):
					c += b[i] * wholeParams[i]
					for j in range(len(b)):
						c += a[i,j] * wholeParams[i] * wholeParams[j]
				c2s.append(c)
			chi2s.append(c2s)
		return chi2s

	def fitParametersForMode(self):
		"""
		Performs the necessary fit dependent on the set mode
		"""
		if self.mode == SMOOTH:
			pass
#			print "Nothing to do for smoothness"
		elif self.mode == AMPL:
			self.fitShapeParameters()
			self.calculateNonShapeParameters()
		elif self.mode == PHASE:
			self.phaseFit()
			self.calculateNonShapeParameters()
		else:
			raise RuntimeError("No mode set")

	def getZeroModeParametersForMode(self):
		"""
		Returns the zero mode parameters obtaiend by the mode
		"""
		if self.mode == SMOOTH:
			params = self.getSmoothnessZeroModeParameters()
		elif self.mode == AMPL or self.mode == PHASE:
			params = self.getZeroModeParameters()
		else:
			raise RuntimeError("No mode set")
		return self.createArtificialSignature(params)


	def evaluateResolvedZeroModeParametersForMode(self, parameters):
		"""
		Will resolve in m and t' bins
		"""
		params = self.reduceArtificialSignatureToOwn(parameters)

		if self.mode == SMOOTH:
			retVal = self.resolvedSmoothnessEvaluation(params)
		elif self.mode == AMPL:
			retVal = self.evaluateZeroModeParameters(params)
		elif self.mode == PHASE:
			retVal =  self.getPhaseChi2s(params, mBinResolved = True)
		else:
			raise RuntimeError("No mode set")
		return retVal

	def setExternalZMP(self, zmp):
		if not len(zmp) == len(self.model):
			raise ValueError("Dimension mismatch (tBin level)")
		for i in range(len(zmp)):
			if not len(zmp[i]) == len(self.model[i]):
				raise ValueError("Dimension mismatch (mBin level)" + str(len(zmp[i])) + " " + str(len(self.model[i])))
		self.externalZMP = zmp
		self.SET("hasExternalZMP")

	def callAmplChi2WithFixedZeroMode(self, params):
		"""
		return a chi2 with given zero-mode parameters
		"""
		if not self.IS("hasExternalZMP"):
			raise Exception("Set the zero-mode parameters first")
		c2 = 0.
		if len(params) == 0:
			raise Exception("Nothing to fit")
		for t,tBin in enumerate(self.model):
			for m,mBin in enumerate(tBin):
				mBin.setShapeParameters(params)
				zmp   = self.externalZMP[t][m]
				a,b,c = mBin.getOwnTheoryABC()
				zerD  = len(zmp)
				newD  = len(b) - zerD
				A     = np.zeros((newD, newD))
				B     = np.zeros((newD))
				C     = c
				for i in range(newD):
					B[i] += b[zerD+i]
					for j in range(newD):
						A[i,j] = a[zerD+i,zerD+j]
				for i in range(zerD):
					C += zmp[i] * b[i]
					for j in range(zerD):
						C += zmp[i]*a[i,j]*zmp[j]
					for j in range(newD):
						B[j] += (a[i,j+zerD] + a[j+zerD,i])*zmp[i]
				cpl = -np.dot(la.pinv(A + np.transpose(A)), B)
				c2  += C + np.dot(B, cpl) + np.dot(cpl, np.dot(A, cpl))
		return c2

	def evaluateZeroModeParametersForMode(self, parameters):
		"""
		This will be returned t'-bin wise
		"""
		params = self.reduceArtificialSignatureToOwn(parameters)

		if self.mode == SMOOTH:
			return self.getSmoothnessChi2s(params)
		elif self.mode == AMPL:
			chi2s = self.evaluateZeroModeParameters(params)
			retVal = []
			for tt in chi2s:
				val = 0.
				for v in tt:
					val += v
				retVal.append(val)
			return retVal
		elif self.mode == PHASE:
			return self.getPhaseChi2s(params)
		else:
			raise RuntimeError("No mode set")

	def getSectorIndex(self, sector):
		"""
		Returns the index of a sectoe, to make it possible to ernter ints of strs as sector specifiers
		"""
		if isinstance(sector, int):
			if sector < 0 or sector >= len(self.sectors):
				raise ValueError("sector '" + str(sector) + "' exceeds range")
			return sector
		for i, s in enumerate(self.sectors):
			if s == sector:
				return i
		raise ValueError("Sector '" + sector + "' invalid")

	def getNDFforMode(self):
		ndfs  = self.model.getNDF()
		ndf   = 0
		nnon  = 0
		first = True
		for tt in ndfs:
			for mm in tt:
				ndf  += mm[0]
				nnon += mm[1] + mm[2]
				if first:
					first = False
					nPar = mm[3]
				else:
					if not nPar == mm[3]:
						raise ValueError("Number of parameters does not match")
		if self.mode == AMPL:
			return ndf - nnon - nPar
		if self.mode == PHASE:
			return ndf/2 - nnon -nPar
		if self.mode == SMOOTH:
			if not nPar == 0:
				raise ValueError("smooth mode has to have no parameters")
			return ndf - nnon
		raise RuntimeError("No mode set")

	def produceResultViewer(self, nonShapeParameters, sector, tBin = -1, noRun = False, plotTheory = False):
		"""
		Produces a result viewer for the given zero mode parameters
		"""
		if tBin == -1:
			if len(self.tBins) == 1:
				tBinString = self.tBins[0]
				tBin = 0
			else:
				raise Exception("tBin not specified")

		for tt in self.model:
			tt.removeZeroModeFromComa()

		self.evaluateZeroModeParametersForMode(self.getZeroModeParametersForMode())
		nSect = self.getSectorIndex(sector)
		lzp = []
		zzz = []
		if plotTheory:
			for i,tt in enumerate(self.model):
				tt.setTheoryFromOwnFunctions(self.nonShapeParameters[i], skipZeroPars = True)

		for ib, mb in enumerate(nonShapeParameters[tBin]):
			for i in range(2*self.model[tBin][ib].nZero):
				lzp.append(nonShapeParameters[tBin][ib][i])
				zzz.append(0.)

		dataInte = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			dataInte.append(histCopy)
		dataReal = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			dataReal.append(histCopy)
		dataImag = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			dataImag.append(histCopy)
		dataPhas = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			dataPhas.append(histCopy)
		corrInte = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			corrInte.append(histCopy)
		corrReal = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			corrReal.append(histCopy)
		corrImag = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			corrImag.append(histCopy)
		corrPhas = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			corrPhas.append(histCopy)
		if plotTheory:
			theoInte = []
			for s in range(len(self.sectors)):
				histCopy = self.histListReal[s].Clone()
				histCopy.Reset()
				theoInte.append(histCopy)
			theoReal = []
			for s in range(len(self.sectors)):
				histCopy = self.histListReal[s].Clone()
				histCopy.Reset()
				theoReal.append(histCopy)
			theoImag = []
			for s in range(len(self.sectors)):
				histCopy = self.histListReal[s].Clone()
				histCopy.Reset()
				theoImag.append(histCopy)
			theoPhas = []
			for s in range(len(self.sectors)):
				histCopy = self.histListReal[s].Clone()
				histCopy.Reset()
				theoPhas.append(histCopy)

		dataCOMA = []
		for s in range(len(self.sectors)):
			histCopy = self.histListReal[s].Clone()
			histCopy.Reset()
			dataCOMA.append(histCopy)

		self.model[tBin].fillHistograms(zzz, dataInte                         )
		self.model[tBin].fillHistograms(zzz, dataReal, mode = REAL            )
		self.model[tBin].fillHistograms(zzz, dataImag, mode = IMAG            )
		self.model[tBin].fillHistograms(zzz, dataPhas, mode = PHASE           )
		self.model[tBin].fillHistograms(zzz, dataCOMA, mode = REIMCORRELATION )

		self.model[tBin].fillHistograms(lzp, corrInte                         )
		self.model[tBin].fillHistograms(lzp, corrReal, mode = REAL            )
		self.model[tBin].fillHistograms(lzp, corrImag, mode = IMAG            )
		self.model[tBin].fillHistograms(lzp, corrPhas, mode = REIMCORRELATION )

		if plotTheory:
			self.model[tBin].fillHistograms(lzp, theoInte, mode = INTENSTHEO)
			self.model[tBin].fillHistograms(lzp, theoReal, mode = REALTHEO)
			self.model[tBin].fillHistograms(lzp, theoImag, mode = IMAGTHEO)
			self.model[tBin].fillHistograms(lzp, theoPhas, mode = PHASETHEO)

		for s in range(len(self.sectors)):
			renormToBinWidth(dataInte[s]    )
			renormToBinWidth(dataReal[s], .5)
			renormToBinWidth(dataImag[s], .5)
			renormToBinWidth(dataCOMA[s],   )
			renormToBinWidth(corrInte[s],   )
			renormToBinWidth(corrReal[s], .5)
			renormToBinWidth(corrImag[s], .5)
			if plotTheory:
				renormToBinWidth(theoInte[s]    )
				renormToBinWidth(theoReal[s], .5)
				renormToBinWidth(theoImag[s], .5)

		if not plotTheory:
			rv = resultViewer([corrInte[nSect], dataInte[nSect]], [corrReal[nSect], dataReal[nSect]], [corrImag[nSect], dataImag[nSect]], [corrPhas[nSect], dataPhas[nSect]], self.startBin, reImCorrel = dataCOMA[nSect], noRun = noRun)
		else:
			rv = resultViewer([corrInte[nSect], dataInte[nSect], theoInte[nSect]], [corrReal[nSect], dataReal[nSect], theoReal[nSect]], [corrImag[nSect], dataImag[nSect], theoImag[nSect]], [corrPhas[nSect], dataPhas[nSect], theoPhas[nSect]], self.startBin, reImCorrel = dataCOMA[nSect], noRun = noRun)

		rv.titleRight = getProperWaveName(self.sectors[nSect])
		rv.tString    = getProperDataSet(self.inFileName, tBinString)
		return rv

def main():
	inFileName       = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_4pp.root"
	sectors          = ["1++0+[pi,pi]1--PiS"]
	tBins            = [0,1]
	startBin         = 20
	stopBin          = 30

	mPi      = 0.13957018

	rhoMass  = ptc.parameter( .77549, "rhoMass" )
	rhoWidth = ptc.parameter( .1491 , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 0, False)

	aa = amplitudeAnalysis(inFileName, sectors, {"1++0+[pi,pi]1--PiS":[rho]}, startBin, stopBin, tBins)
	aa.loadData()
	aa.finishModelSetup()
	aa.fitShapeParameters()
	aa.calculateNonShapeParameters()
	params = aa.getZeroModeParameters()
	c2s = aa.getPhaseChi2s(params)
	print c2s
	for i in range(len(params)):
		if not i ==0:
			continue
		for j in range(len(params[i])):
			for k in range(len(params[i][j])):
				params[i][j][k] = 1.345
	c2s = aa.getPhaseChi2s(params)
	print c2s
	vl = 0.
	for v in c2s:
		vl += v
#	print "Total:",vl
#	sys.exit(0)

	parSmooth = aa.getSmoothnessZeroModeParameters()
	aa.fitShapeParameters()
	aa.calculateNonShapeParameters()
	params = aa.getZeroModeParameters()

	fif = sumUp(aa.evaluateZeroModeParameters(params))
	sif = sumUp(aa.evaluateZeroModeParameters(parSmooth))

	fis = sumUp(aa.getSmoothnessChi2s(params))
	sis = sumUp(aa.getSmoothnessChi2s(parSmooth))

	print "fit in fit"
	print fif
	print
	print "smooth in fit"
	print sif
	print (sif-fif)/fif
	print
	print "fit in smooth"
	print fis
	print (fis-sis)/sis
	print
	print "smooth in smooth"
	print sis



	sys.exit(0)

	aaR = amplitudeAnalysis(inFileName, sectors, {"1++0+[pi,pi]1--PiS":[rho]}, startBin, stopBin, tBins, sectorRangeMap = {'1++0+[pi,pi]1--PiS':(.65,1.)})
	aaR.loadData()
	aaR.finishModelSetup()
#	aaR.phaseFit()
	aaR.fitShapeParameters()
	aaR.calculateNonShapeParameters()
	paramsR = aaR.getZeroModeParameters()
	chi2s = aa.evaluateZeroModeParameters(params)
	chi2r = aa.evaluateZeroModeParameters(paramsR)

	chi2rs = aaR.evaluateZeroModeParameters(paramsR)
	chi2rr = aaR.evaluateZeroModeParameters(params)


	print sumUp(chi2s)  ,  sumUp(chi2r)
	print sumUp(chi2rs) , sumUp(chi2rr)




	count =  0
	totalChi2 = 0.
	with open("resultsFile", 'w') as  out:
		for i in range(len(chi2s)):
			for j in range(len(chi2s[0])):
				totalChi2 += chi2s[i][j]
				print (chi2s[i][j] - chi2r[i][j])#/chi2s[i][j]
				out.write(str(count) +  ' ' + str((chi2s[i][j] - chi2r[i][j])) + ' ' + str((chi2s[i][j] - chi2r[i][j])/chi2s[i][j]) + '\n')
				count += 1
	print totalChi2, aa.chi2
#	chi2 = 0.

#	for tt in chi2s:
#		for mm in tt:
#			chi2 += mm
#	print chi2, aa.chi2
if __name__ == "__main__":
	main()
