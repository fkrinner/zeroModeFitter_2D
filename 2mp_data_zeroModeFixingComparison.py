import matplotlib
matplotlib.use("Agg")

from analysisClass import amplitudeAnalysis
import parameterTrackingParameterizations as ptc
import parameterizationClasses as pc
from fixedparameterizationPaths import getFileNameForSector
from modes import PHASE, AMPL, SMOOTH, NONE
from utils import sumUp, weightedSum, cloneZeros, checkLaTeX
import sys
import ROOT
import numpy as np
from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0,m1500,G1500, mF2prime, GF2prime
from rootfabi import root_open
import LaTeX_strings

import consistencyUtils as cu
from studyFileNames import fileNameMap

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from random import randint, random
import scipy

from LaTeX_strings import unCorrected_string, weightedAVG_string

def doFunctionFit(inFileName, funcs, startBin, stopBin, tBins, sectorRangeMap, referenceWave = "", writeResultToFile = None, acv = None, zeroModeParameters = None, ifn = None):
	sector  = "2-+0+[pi,pi]2++PiS"
	sectors = ["2-+0+[pi,pi]0++PiD","2-+0+[pi,pi]1--PiP","2-+0+[pi,pi]1--PiF","2-+0+[pi,pi]2++PiS"]

	print "La wave du reference e >>>",referenceWave,"<<<",sectorRangeMap
	fitter = amplitudeAnalysis(inFileName, sectors, {sector:funcs}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitter.loadData(loadIntegrals = True, referenceWave = referenceWave)
	fitter.finishModelSetup()
	fitter.mode = AMPL

	parameters = []
	for f in funcs:
		for p in f.returnParameters():
			parameters.append(p)

	fitter.initMinuitFunction(parameters)

	fitter.removeZeroModeFromComa()
	fitter.removeGlobalPhaseFromComa()


	if acv is not None:
		print "adding ACV of",acv
		fitter.addComaValueForZeroMode(acv)

	if ifn is None:
		if zeroModeParameters is None:
			fitter.fitShapeParameters()
			fitter.calculateNonShapeParameters()
			chi2 = fitter.chi2
			ndf  = fitter.getNDFforMode()
		else:
			x,err,chi2,ndf = fitter.fitShapeParametersForBinRange([p.value for p in parameters], [0], range(stopBin-startBin), zeroModeParameters = zeroModeParameters)
			fitter.calculateNonShapeParametersForZeroModeParameters(zeroModeParameters)
			for i,p in enumerate(parameters):
				print p
		if writeResultToFile is not None:
			with open(writeResultToFile, 'w') as outFile:
				outFile.write("- - - - parameters - - - -\n")
				for p in parameters:
					outFile.write(str(p)+'\n')
				outFile.write(" - - - - - fit - - - - -\nchi2/NDF: "+str(chi2)+"/"+str(ndf)+"="+str(chi2/ndf))
	else:
		params = []
		errs   = []
		with open(ifn, 'r' ) as inFile:
			inPars = False
			for line in inFile.readlines():
				if "- - - - - fit - - - - -" in line:
					break
				elif inPars:
					params.append(float(line.split()[2]))
					errs.append(float(line.split()[4]))
				elif "- - - - parameters - - - -" in line:
					inPars = True
		if not len(params) == len(parameters):
			raise IndexError("Parameter size mismatch")
		for i in range(len(params)):
			parameters[i].value = params[i]
			parameters[i].error = errs[i]
		if zeroModeParameters is None:
			raise ValueError
		else:
			fitter.setZeroModeParameters(zeroModeParameters)
		fitter.model.setBinsToEvalueate([0], range(stopBin-startBin))
		fitter.chi2 = fitter.model.fixedZMPchi2(params)
		fitter.fitParameters = params
		fitter.SET('hasFitResult')
#		fitter.calculateNonShapeParametersForZeroModeParameters(zeroModeParameters)
		fitter.calculateNonShapeParameters()
		
	return fitter

def doF0Fit(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", deg = 0, writeResultToFile = None):
	fixedAMPD_sub = pc.fixedParameterization("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_2mp0pSigmaPiD")
	f0Mass   = ptc.parameter(mF0,"f0Mass")
	f0g1     = ptc.parameter(g1, "g1")
	f0g2     = ptc.parameter(g2, "g2")
#	f0g1.lock = True
#	f0g2.lock = True
	f0       = ptc.flatte([f0Mass, f0g1, f0g2], mPi, mPi, mPi, mK, mK, 0, 1, False)
	polyPars = []
	for i in range(deg):
		polyPars.append(ptc.parameter(0.1, "reC"+str(i+1)))
		polyPars.append(ptc.parameter(0.1, "imC"+str(i+1)))
	poly  = ptc.complexPolynomial(deg, polyPars)
	fitF0 = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]0++PiD":[f0,fixedAMPD_sub]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitF0.loadData(referenceWave = referenceWave)
	fitF0.finishModelSetup()
	fitF0.fitShapeParameters()
	fitF0.calculateNonShapeParameters()
	fitF0.mode = AMPL
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(f0Mass.value) + ' ' + str(f0Mass.error) + ' ' + str(f0g1.value) + ' ' + str(f0g1.error) + ' ' + str(f0g2.value) + ' ' + str(f0g2.error) 
			for i in range(deg):
				rasultString += ' ' + str(polyPars[i].value) + ' ' + str(polyPars[i].error)
			resultString +=  "\n"
			outFile.write(resultString)
	return fitF0

def doFit1500(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"2-+0+[pi,pi]0++PiD": (1.15, 1.75)}, referenceWave = "", writeResultToFile = None):
	fixedAMPD_sub = pc.fixedParameterization("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_2mp0pSigmaPiD")
#	fixedAMPD_sub = pc.constant()
	f0mass  = ptc.parameter( m1500,  "1500mass" )
	f0width = ptc.parameter( G1500 , "1500width")
	f0_1500 = ptc.relativisticBreitWigner([f0mass,f0width], mPi, mPi, mPi, 0, 0, False)
	fit1500 = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]0++PiD":[f0_1500, fixedAMPD_sub]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fit1500.loadData(referenceWave = referenceWave)
	fit1500.finishModelSetup()
	fit1500.fitShapeParameters()
	fit1500.calculateNonShapeParameters()
	fit1500.mode = AMPL
#	fit1500.removeGlobalPhaseFromComa()
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(f0mass.value) + ' ' + str(f0mass.error) + ' ' + str(f0width.value) + ' ' + str(f0width.error) + "\n"
			outFile.write(resultString)
	return fit1500

def doFitRhoP(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(rhoMass.value) + ' ' + str(rhoMass.error) + ' ' + str(rhoWidth.value) + ' ' + str(rhoWidth.error) + "\n"
			outFile.write(resultString)
	return fitRho

def doFitRhoF(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiF":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(rhoMass.value) + ' ' + str(rhoMass.error) + ' ' + str(rhoWidth.value) + ' ' + str(rhoWidth.error) + "\n"
			outFile.write(resultString)
	return fitRho

def doFitBothRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rhoP = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	rhoF = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiP":[rhoP], "2-+0+[pi,pi]1--PiF" : [rhoF]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(rhoMass.value) + ' ' + str(rhoMass.error) + ' ' + str(rhoWidth.value) + ' ' + str(rhoWidth.error) + "\n"
			outFile.write(resultString)
	return fitRho

def doFitF2(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None):
	f2Mass  = ptc.parameter( mF2,  "f2Mass" )
	f2Width = ptc.parameter( GF2 , "f2Width")
	f2 = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, False)
	fitF2 = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]2++PiS":[f2]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitF2.loadData(referenceWave = referenceWave)
	fitF2.finishModelSetup()
	fitF2.fitShapeParameters()
	fitF2.calculateNonShapeParameters()
	fitF2.mode = AMPL
#	fitF2.removeGlobalPhaseFromComa()
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(f2Mass.value) + ' ' + str(f2Mass.error) + ' ' + str(f2Width.value) + ' ' + str(f2Width.error) + "\n"
			outFile.write(resultString)
	return fitF2

f2PrimeFit        = None # Ugly but convenient
def fitF2prime(fitBin, outFileName, inFileName, startBin, stopBin, tBins, pars, zeroModeParameters, plotNameBase = "", referenceWave = ""):
	global f2PrimeFit
	if not f2PrimeFit:
		f2Mass  = ptc.parameter( pars[0],  "f2Mass")
		f2Width = ptc.parameter( pars[1], "f2Width")
		f2Mass.lock  = True
		f2Width.lock = True
		f2 = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, False)
		primeMass   = ptc.parameter( mF2prime, "primeMass" )
		primeWidth  = ptc.parameter( GF2prime, "primeWidth")
		prime       = ptc.relativisticBreitWigner([primeMass,primeWidth], mPi, mPi, mPi, 2, 0, False)
		f2PrimeFit = amplitudeAnalysis(inFileName, ["2-+0+[pi,pi]2++PiS"], {"2-+0+[pi,pi]2++PiS":[f2,prime]}, startBin, stopBin, tBins)
		f2PrimeFit.loadData(referenceWave = referenceWave)
		f2PrimeFit.finishModelSetup()
		f2PrimeFit.fitShapeParameters()
		f2PrimeFit.calculateNonShapeParameters()
		f2PrimeFit.mode = AMPL
	binIndex = fitBin+startBin
	with open(outFileName, 'a') as outFile:
		outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
		startValueOffset = 0.01
		exceptCount      = 0
		while True:
#			try:
			if True:
				f2PrimeFit.initMinuitFunction(f2PrimeFit.getParameters()[2:])
				x,err,c2,ndf = f2PrimeFit.fitShapeParametersForBinRange([mF2prime+startValueOffset,GF2prime+startValueOffset], [0],[fitBin], zeroModeParameters = zeroModeParameters)
				break
#			except:
#				print "Fitter exception encountered"
#				startValueOffset += 0.001
#				exceptCount      += 1
#				if exceptCount > 3:
#					print "Too many failed attempts in bin "+str(fitBin)+": "+str(exceptCount)
#					raise Exception
#					x, err = [0.,0.],[0.,0.]
#					break

		outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
		outFile.write(' ' + str(c2/ndf) + '\n')
		if not plotNameBase == "":
			f2PrimeFit.calculateNonShapeParameters()
			f2primeRV = f2PrimeFit.produceResultViewer(zeroModeParameters,"2-+0+[pi,pi]2++PiS", noRun = True, plotTheory = True)
			f2primeRV.plotData = True
			f2primeRV.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])

def doF0phase(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"2-+0+[pi,pi]0++PiD" : (0.34, 2*mK - .04)}, referenceWave = ""):
	pipiSfileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_2mp0pSigmaPiD"
	pipiSw = pc.fixedParameterization(pipiSfileName, polynomialDegree  = 0, complexPolynomial = False)
	waveModel = {"2-+0+[pi,pi]0++PiD": [pipiSw]}
	fitPiPiSshape = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitPiPiSshape.loadData(referenceWave = referenceWave)
	fitPiPiSshape.finishModelSetup()
	fitPiPiSshape.phaseFit()
	fitPiPiSshape.mode = PHASE
#	fitPiPiSshape.removeGlobalPhaseFromComa()
	fitPiPiSshape.unifyComa()
	return fitPiPiSshape

alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	waveModel = {}
	for n,sector in enumerate(sectors):
		model = []
		fileNames = getFileNameForSector(sector, False, False)
		print fileNames
		for fn in fileNames:
			param = pc.fixedParameterization(fn, polynomialDegree  = 0, complexPolynomial = False)
			model.append(param)
		waveModel[sector] = model
	fixedShapes = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fixedShapes.loadData(loadIntegrals = True, referenceWave = referenceWave)
	fixedShapes.finishModelSetup()
	fixedShapes.fitShapeParameters()
	fixedShapes.calculateNonShapeParameters()
	fixedShapes.mode = AMPL
#	fixedShapes.removeGlobalPhaseFromComa()
	return fixedShapes

def doSmooth(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	waveModel = {}
	smooth = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	smooth.loadData(referenceWave = referenceWave)
	smooth.finishModelSetup()
	smooth.mode = SMOOTH
#	smooth.removeGlobalPhaseFromComa()
#	smooth.removeZeroModeFromComa()
	return smooth

def weightedParametersSum(evl, selfEvl, params):
	"""
	Calculates the weight of each method as: weight[n] = 1./sum{(evl[m,n] - selfEvl[m])/selfEvl[m] , m != n}
	"""
	weights     = {}
	totalWeight = 0.
	for n in selfEvl:
		diffSum = 0.
		for m in selfEvl:
			if m == n:
				continue
			diffSum += (evl[m,n] - selfEvl[m])/selfEvl[m]
		weight = 1./diffSum
		totalWeight += weight
		weights[n] = weight
	for m in weights:
		weights[m] /= totalWeight
	return weightedSum(weights, params)

def main():
	checkLaTeX()
	style = modernplotting.mpplot.PlotterStyle()
#	style.p2dColorMap = 'ocean_r'
#	style.p2dColorMap = 'YlOrRd'
	style.p2dColorMap = 'Reds'
#	referenceWave     = ""	

	sectors          = ["2-+0+[pi,pi]0++PiD", "2-+0+[pi,pi]1--PiP", "2-+0+[pi,pi]1--PiF", "2-+0+[pi,pi]2++PiS"]
	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))
	if len(sys.argv) > 2 and False:
		study = sys.argv[2]
		studyAdder = "_"+study
	else:
		study      = "std11"
		studyAdder = ""
	print "Study: "+study

	inFileName = fileNameMap[study]

	tBins          = [tBin]
	rhoRange       = None
	f2Range        = None
	sectorRangeMap = {}

	if rhoRange is None:
		rhoRangeString = ""
	else:
		sectorRangeMap["2-+0+[pi,pi]1--PiP"] = (0.,rhoRange)
		sectorRangeMap["2-+0+[pi,pi]1--PiF"] = (0.,rhoRange)
		rhoRangeString = "_range"+str(rhoRange)
	if f2Range is None:
		f2RangeString = ""
	else:
		sectorRangeMap["2-+0+[pi,pi]2++PiS"] = (0., f2Range)
		f2RangeString = "_range"+str(f2Range)

	startBin         = 11
	stopBin          = 50

#	startBin         = 40
#	stopBin          = 50

	ifn              = None
	for a in sys.argv:
		if a.startswith("mBins"):
			startBin = int(a.split("-")[1])
			stopBin  = int(a.split("-")[2])
		if a.startswith("ifn="):
			ifn = a[4:]

	methodBinRanges = {
	                   "fitF2"        : (22, 50),
	                   "fitRhoF"      : (22, 50),
	                   "fixedShapeF2" : (22, 50)}
#	methodBinRanges = {} # Override here 



	allMethods       = {}
	methodStrings    = {}
	shortlabels      = {  "fixedShapeF0"      : r"$\text{fix}_{f_0}^{~}$",
	                      "fixedShapeRhoP"    : r"$\text{fix}_\rho^{P}$",
	                      "fixedShapeRhoF"    : r"$\text{fix}_\rho^{F}$",
	                      "fixedShapeBothRho" : r"$\text{fix}_\rho^{2}$",
	                      "fixedShapeF2"      : r"$\text{fix}_{f_2}$",
	                      "fixedShapes"       : r"$\text{fix}_\text{all}^{~}$",
	                      "pipiS"             : r"$\phi_{[\pi\pi]_S}^{~}$",
	                      "fitRhoP"           : r"$\text{fit}_\rho^{P}$",
	                      "fitRhoF"           : r"$\text{fit}_\rho^{F}$",
	                      "fitBothRho"        : r"$\text{fit}_\rho^{2}$",
	                      "fitF2"             : r"$\text{fit}_{f_2}$",
	                      "smooth"            : r"smooth"}


#       # - - - - --- Start here with the model builting --- - - - - #       #

	nPol = 3
	nCmx = 0
	for a in sys.argv:
		if a.startswith("nPol"):
			nPol = int(a[4:])
		if a.startswith("nCmx"):
			nCmx = int(a[4:])
	f2Re0     = mF2**2
	f2Im0     = mF2*GF2

	mPrime = 1.55
	Gprime =  .2

	mPPrime = 1.91
	Gpprime = .3


	f2Mass    = ptc.parameter(mF2, "f2_mass")
	f2Width   = ptc.parameter(GF2, "f2_width")
	f2        = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, False)

	poleReal   = ptc.parameter(f2Re0, "f2Re")
	poleImag   = ptc.parameter(f2Im0, "f2Im")

	pPoleReal = ptc.parameter(mPrime**2,      "f2pRe")
	pPoleImag = ptc.parameter(mPrime*Gprime, "f2pIm")

	ppPoleReal = ptc.parameter(mPPrime**2, "f2ppRe")
	ppPoleImag = ptc.parameter(mPPrime*Gpprime, "f2ppIm")


#	poleRealPrime = ptc.parameter(mPrime**2, "rhoRePrime")
#	poleImagPrime = ptc.parameter(mPrime*Gprime, "rhoImPrime")

	seedint = randint(0,10000)

	polyDeg_po = 3
	for a in sys.argv:
		if a.startswith("pdpo"):
			polyDeg_po = int(a[4:])
	
	params     = [poleReal,poleImag,pPoleReal,pPoleImag,ppPoleReal,ppPoleImag]
	params     = params[:2*nPol]
	for i in range(nCmx):
		params.append(ptc.parameter(2.,"comPolRe_"+str(i)))
		params.append(ptc.parameter(0.,"comPolIm_"+str(i)))
		params.append(ptc.parameter(1.,"comPolG_"+str(i)))

	for d in range(polyDeg_po):
		params.append(ptc.parameter(2*random()-1., "c_"+str(d)))
	Kmatrix    = ptc.simpleOneChannelKmatrix(params, nPol, polyDeg_po, 4*mPi**2,nComplex = nCmx)
	useCM      = False
	for a in sys.argv:
		if a == "CM":
			useCM = True
		if a == "rho":
			useCM = False

	Kmatrix.use_CM = useCM

	pPolyDeg3  = 2
	pPolyDeg2  = 5
	for a in sys.argv:
		if a.startswith("m3pol"):
			pPolyDeg3 = int(a[5:])
		if a.startswith("m2pol"):
			pPolyDeg2 = int(a[5:])


	params     = []
	for d in range(pPolyDeg2):
		for e in range(pPolyDeg3+1):
			params.append(ptc.parameter(2*random()-1., "c_"+str(d+1)+"_"+str(e)))
	pPoly      = ptc.twoDimensionalRealPolynomial(pPolyDeg2, pPolyDeg3, params, baseExponent = 2) # baseExponent = 2: polynomial in s
	func       = ptc.multiply([Kmatrix, pPoly])	

	model      = [func]
#	model      = [rho]

	acv = None

#       # - - - - --- Stop the model building here --- - - - - #       #
	zeroModeParameters = None
	fixZeroMode        = True
#	print sectorRangeMap
#	return 
	if fixZeroMode:
		fixSectors = ["2-+0+[pi,pi]1--PiP", "2-+0+[pi,pi]1--PiF", "2-+0+[pi,pi]2++PiS"]
		fixRangeMap = {
			"2-+0+[pi,pi]1--PiP" : (0.,1.12), 
			"2-+0+[pi,pi]1--PiF" : (0.,1.12), 
			"2-+0+[pi,pi]2++PiS" : (0.,1.27) # Use only left side of f_2(1270)
		}

		fixedShapes        = doFixedShapes(inFileName, fixSectors, startBin, stopBin, tBins, referenceWave = referenceWave, sectorRangeMap = fixRangeMap)
		zeroModeParameters = fixedShapes.getZeroModeParametersForMode()

#		RV = fixedShapes.produceResultViewer(zeroModeParameters,"1-+1+[pi,pi]1--PiP", noRun = True, plotTheory = True)
#		RV.plotData = True
#		for b in range(startBin, stopBin):
#			plotNameBase = "./Kmatrix_plots/1mp1p1mmPiP_<mode>_"+str(b)+"_"+str(tBin)+".pdf"
#			RV.writeBinToPdf(b, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
#		return
	if useCM:
		ps = "CM"
	else:
		ps = "rho"
	if rhoRange is None:
		rrs = ""
	else:
		rrs = "range"+str(rhoRange)+'_'
	if nPol == 1:
		nps = ""
	else:
		nps = "nPol"+str(nPol)+"_"
	if nCmx > 0:
		print "Using complex pole"
		nps += "nCmx"+str(nCmx)+"_"

	resultFile  =  "./complexF2poles/2mpF2_Kmatrix_"+rrs+nps+ps+"_kPol"+str(polyDeg_po-1)+"_pPol"+str(pPolyDeg3)+"-"+str(pPolyDeg2)+"_t"+str(tBin)+"_m"+str(startBin)+'-'+str(stopBin)+'_'+str(seedint)+".dat"
	fitter      = doFunctionFit(inFileName, model, startBin, stopBin, tBins, sectorRangeMap, referenceWave = referenceWave, acv = acv, zeroModeParameters = zeroModeParameters, writeResultToFile = resultFile, ifn = ifn)

#       # - - - - --- Start the evaluations here --- - - - - #       #
	nBinsPlot = 1000
	def fPlot(v):
		return v.imag
#	hist = ROOT.TH2D("hhh","hhh", nBinsPlot, -.25, 6.25,nBinsPlot, -1., 1.)
#	for iX in range(nBinsPlot):
#		x = hist.GetXaxis().GetBinCenter(iX+1)
#		for iY in range(nBinsPlot):
#			y = hist.GetYaxis().GetBinCenter(iY+1)			
#			s = x+1.j*y
#			val = Kmatrix.complexCall(s)
#			hist.SetBinContent(iX+1, iY+1,fPlot(val))
#	hist.Draw("COLZ")
#	raw_input("press <enter> to go to the secont sheet")
	Kmatrix.secondSheet = True
#	for iX in range(nBinsPlot):
#		x = hist.GetXaxis().GetBinCenter(iX+1)
#		for iY in range(nBinsPlot):
#			y = hist.GetYaxis().GetBinCenter(iY+1)			
#			s = x+1.j*y
#			val = Kmatrix.complexCall(s)
#			hist.SetBinContent(iX+1, iY+1, fPlot(val))
#	hist.Draw("COLZ")
#	res = scipy.optimize.minimize(Kmatrix.absInverse,[f2Re0,f2Im0])
#	print res.x,"pole position"
#	mfv = res.fun
#	resSting = str(res.fun)+ " function value should be zero"
#	print resSting
#	BWstring = "BW par: "+str(abs(res.x[0])**.5)+" "+str(abs(res.x[1])/abs(res.x[0])**.5)+" (all absolute values)"
#	print BWstring

#	if ifn is None:
#		print "= = = = = = = = = Starting BW error ersimation = = = = = = = = = "
#		nPoints     = 1000
#		poleMean    = res.x
#		poleSamples = []
#		i           = 0
#		failCount   = 0
#		while i < nPoints:
#			pts = np.random.multivariate_normal(fitter.fitParameters, fitter.MINUITcoma)
#			fitter.MINUIT_function(pts) # Call the function once to set parameters inside
#	#		fitter.model[0].setParametersAndErrors(pts, fitter.MINUITerrs)
#			res = scipy.optimize.minimize(Kmatrix.absInverse,poleMean)
#			if abs(res.fun) > 100*mfv:
#				print "No more pole found (mfv = "+str(mfv)+") : fval = "+str(res.fun)
#				failCount += 1
#				if failCount > nPoints:
#					print "Failed to find poles too often.... abort"
#					return
#				continue
#	#			raise ValueError("No more pole found: fval = "+str(res.fun))
#
#			poleSamples.append(res.x)
#			i+=1
#	#		print i,"Marker to find the PRINT 57473M3N7"
#		meanPole = [0.,0.]
#		for p in poleSamples:
#			meanPole[0] += p[0]
#			meanPole[1] += p[1]
#		meanPole[0] /= len(poleSamples)
#		meanPole[1] /= len(poleSamples)
#		poleComa = [[0.,0.],[0.,0.]]
#		for p in poleSamples:
#			poleComa[0][0] += (p[0]-meanPole[0])**2
#			poleComa[0][1] += (p[0]-meanPole[0])*(p[1]-meanPole[1])
#			poleComa[1][0] += (p[1]-meanPole[1])*(p[0]-meanPole[0])
#			poleComa[1][1] += (p[1]-meanPole[1])**2
#		poleComa[0][0] /= len(poleSamples)-1
#		poleComa[0][1] /= len(poleSamples)-1
#		poleComa[1][0] /= len(poleSamples)-1
#		poleComa[1][1] /= len(poleSamples)-1
#		comaString      = str(poleComa)
#		print " - - - - - - le compaire pramaitre  - - - - - - "
#		print meanPole, poleMean
#		print " - - - - - - le compaire pramaitre  - - - - - - "
#		print poleComa
#		print "= = = = = = = = = Finished BW error ersimation = = = = = = = = = "
#		mF2P = 1.9
#		GF2P =  .277
#		res  = scipy.optimize.minimize(Kmatrix.absInverse,[mF2P**2,mF2P*GF2P])
#		print res.x,"pole position"
#		resSting = str(res.fun)+ " function value should be zero"
#		print resSting
#		BWstring = "BW' par: "+str(abs(res.x[0])**.5)+" "+str(abs(res.x[1])/abs(res.x[0])**.5)+" (all absolute values)"
#		print BWstring
#		with open(resultFile,'a') as outFile:
#			outFile.write('\n'+BWstring+" "+resSting+"\ncoma "+comaString)
#	#		outFile.write('\n'+BWfitString)
#	#       # - - - - --- Start the evaluations here --- - - - - #       #
#		doPlots = False
#		for a in sys.argv:
#			if a == "plot":
#				doPlots = True
#	if ifn is not None:
#		doPlots = True # Set plots by default, if an inFile is given
#
#	if doPlots:
#		RV = fitter.produceResultViewer(zeroModeParameters,"2-+0+[pi,pi]2++PiS", noRun = True, plotTheory = True)
#		RV.plotData = True
#		for b in range(startBin, stopBin):
#			plotNameBase = "./Kmatrix_plots/2mp0p2ppPiS_<mode>_"+str(b)+"_"+str(tBin)+"_"+str(seedint)+".pdf"
#			RV.writeBinToPdf(b, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
#	raw_input("press <enter> to exit")
#	return
##################################################################################
#################################################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################


#	fit1500 = doFit1500(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "2mp_f0_1500_massesAndWidths_global.dat")

	print "Start with fixed shape f0"
	fixedShapeF0 = doFixedShapes(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"


#	fixedShapeAllButF0 = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, referenceWave = referenceWave)
#	hists = fixedShapeAllButF0.makeTheoryTotals()
#	with root_open("theo_hists_2mp0p_t"+str(tBin)+".root", "RECREATE"):
#		for tb in hists:
#			for h in tb:
#				h.Write()
#	return

#	print "Start with fixed shape rhoP"
#	fixedShapeRhoP = doFixedShapes(inFileName, [sectors[1]], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeRhoP"] = fixedShapeRhoP
#	print "Finished with fixed shape rhoP"

#	print "Start with fixed shape rhoF"
#	fixedShapeRhoF = doFixedShapes(inFileName, [sectors[2]], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeRhoF"] = fixedShapeRhoF
#	print "Finished with fixed shape rhoF"

#	print "Start with fixed shape bothRho"
#	fixedShapeBothRho = doFixedShapes(inFileName, sectors[1:3], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeBothRho"] = fixedShapeBothRho
#	print "Finished with fixed shape bothRho"

#	print "Start with fixed shape f2"
#	fixedShapeF2 = doFixedShapes(inFileName, [sectors[3]], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeF2"] = fixedShapeF2
#	print "Finished with fixed shape f2"

	print "Start with fixed shapes"
	fixedShapes = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"
	fullSig = fixedShapes.getZeroModeSignature()

#	totalHists = fixedShapes.getTotalHists(cloneZeros(fixedShapes.getZeroModeParametersForMode()))
#	with root_open("./totals_2mp_noCorr.root", "UPDATE") as out:
#		for t in totalHists:
#			for m in t:
#				m.Write()
#	return

##	print "Start with phase"
##	fitPiPiSshape = doF0phase(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave)
##	allMethods["pipiS"] = fitPiPiSshape
##	print "Finished with phase"

	fitRhoP = doFitRhoP(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "rhoMassesAndWidths_2-+0+1--P_global"+rhoRangeString+".dat", sectorRangeMap = sectorRangeMap)
	fitRhoF = doFitRhoF(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "rhoMassesAndWidths_2-+0+1--F_global"+rhoRangeString+".dat", sectorRangeMap = sectorRangeMap)

	print "Start with fitting bothRho"
	fitBothRho = doFitBothRho(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "rhoMassesAndWidths_2-+0+1--_global"+rhoRangeString+".dat")
	allMethods["fitBothRho"] = fitBothRho
	print "Finished with fitting bothRho"

	print "Start with fitting f2"
	fitF2 = doFitF2(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, referenceWave = referenceWave, writeResultToFile = "f2MassesAndWidths_2-+0+2++_global"+f2RangeString+".dat")
	allMethods["fitF2"] = fitF2
	print "Finished with fitting f2"

	if stopBin - startBin > 1:
		print "Start with smooth"
		smooth = doSmooth(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
		allMethods["smooth"] = smooth
		print "Finished with smooth"

	if "fixedShapeRhoF" in allMethods:
		allMethods["fixedShapeRhoF"].setZeroModeSignature(fullSig,1)
	if "fitRhoF" in allMethods:
		allMethods["fitRhoF"].setZeroModeSignature(fullSig,1)

	diffsFull, resolvedWA, nonResolvedWA, comps, resDiffs, nonResDiffs, resolvedDiffsFull,noCorrDiffs = cu.doAllComparisons(allMethods, startBin, methodBinRanges)
#	print resolvedDiffsFull

	RV = fixedShapes.produceResultViewer(resolvedWA,"2-+0+[pi,pi]2++PiS", noRun = True, plotTheory = True)
	RV.plotData = True
	RV.writeBinToPdf(startBin, stdCmd = ["2mp0p2ppS_2D_"+str(tBin)+".pdf","", [], "", []])
	return

	from math import isnan
	for pair in resolvedDiffsFull:
		with  modernplotting.toolkit.PdfWriter('./resolvedDiffPlots/2mp_'+pair[0]+"_"+pair[1]+"_"+str(tBin)+studyAdder+".pdf") as pdfOutput:
			plot  = style.getPlot1D()
			line  = [0.000000001]*len(resolvedDiffsFull[pair][0])
			line2 = [0.000000001]*len(resolvedDiffsFull[pair][0])
			one   = [1.]*len(resolvedDiffsFull[pair][0])
			xAxis = [ .5 + 0.04*(startBin + i) for i in range(len(resolvedDiffsFull[pair][0]))]
			for i,v in enumerate(resolvedDiffsFull[pair][0]):
				if not isnan(v) and not v <= 0.:
					line[i] = v
				else:
					line[i] = 0.000000001

			if not pair[1] == "WAres" and not pair[1] == "WAnon":
				for i,v in enumerate(resolvedDiffsFull[pair[1],pair[0]][0]):
					if not isnan(v) and not v <= 0.:
						line2[i] = v
					else:
						line2[i] = 0.000000001

			plot.setYlog()
			plot.plot(xAxis, line)
			plot.plot(xAxis, one)
			plot.plot(xAxis, line2)
			plot.setYlim(0.00001, 10000)
			pdfOutput.savefigAndClose()

	studyList = []
	for m in allMethods:
		studyList.append(m)
	studyList.sort()

	style.titleRight = r"$2^{-+}0^+$"
	style.titleLeft  = LaTeX_strings.tBins[tBin]

#	fixedShapesRhosAndF2 = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, referenceWave = referenceWave)

#	with  modernplotting.toolkit.PdfWriter("compositions_2mp_"+str(tBin)+studyAdder+".pdf") as pdfOutput:
#		plot = style.getPlot1D()
#		for m in studyList:
#			line  = [0.]*len(comps[m][0])
#			xAxis = [ .5 + 0.04*(startBin + i) for i in range(len(comps[m][0]))]
#			break
#		count = 0
#		for m in studyList:
#			newLine = line[:]
#			for i in range(len(comps[m][0])):
#				newLine[i] += comps[m][0][i]
#			plot.axes.fill_between(xAxis, line, newLine, facecolor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, 0.1*count))
#			count += 1
#			line = newLine
#		plot.setYlim(0.,1.)
#		plot.setXlim(xAxis[0], xAxis[-1])
#		pdfOutput.savefigAndClose()

	hist = ROOT.TH2D("hist","hist", len(studyList)+2, 0, len(studyList)+2, len(studyList), 0, len(studyList))

	for i,m in enumerate(studyList):
		for j,n in enumerate(studyList):
			hist.SetBinContent(i+1, j+1, diffsFull[n,m])
	for i,m in enumerate(studyList):
		hist.SetBinContent(len(studyList)+1, i+1, noCorrDiffs[m])
		hist.SetBinContent(len(studyList)+2, i+1, resDiffs[m])

	axolotl = []
	for i,study in enumerate(studyList):
		axolotl.append(shortlabels[study])
#		axolotl.append(alphabet[i])

#	with modernplotting.toolkit.PdfWriter("studies_2mp_"+str(tBin)+studyAdder+".pdf") as pdfOutput:
#		plot = style.getPlot2D()
#		plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(len(studyList)+2)])
#		plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(len(studyList))])
#		studyPlotter.makeValuePlot(plot, hist)
		
#		plot.axes.set_yticklabels(axolotl)
#		axolotl.append(unCorrected_string)
#		axolotl.append(weightedAVG_string)
#		plot.axes.set_xticklabels(axolotl, rotation = 90)
#		plot.setZlim((0.,1.))

#		pdfOutput.savefigAndClose()

#	with open("studies_2mp_"+str(tBin)+studyAdder+".txt",'w') as out:
#		for axl in axolotl:
#			out.write(axl + ' ')
#		out.write("\n")
#		for i in range(hist.GetNbinsX()):
#			for j in range(hist.GetNbinsY()):
#				out.write(str(hist.GetBinContent(i+1, j+1)) + ' ')
#			out.write('\n')

	doF0fitGlobal = False
	if doF0fitGlobal:
		deg = 0
#		f0fit = doF0Fit(inFileName, sectors[:1], startBin, stopBin, tBins, sectorRangeMap = {"2-+0+[pi,pi]0++PiD":(0.9,1.1)}, referenceWave = referenceWave, deg = deg)
		f0fit = doF0Fit(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave, deg = deg)
		startPars = [mF0,g1,g2]
		for _ in range(deg):
			startPars.append(0.)
			startPars.append(0.)
		x,err,c2,ndf = f0fit.fitShapeParametersForBinRange(startPars,[0],range(stopBin-startBin), zeroModeParameters = resolvedWA)
		f0fit.setShapeParameters(x,err,resolvedWA)
#		x,err,c2,ndf = f0fit.fitShapeParametersForBinRange([mF0],[0],range(stopBin-startBin), zeroModeParameters = resolvedWeightedSum)
#		print x
		s = 0 # only one sector??
		rv = f0fit.produceResultViewer(resolvedWA,s, noRun = True, plotTheory = True)
		for b in range(startBin, stopBin):
			intensName = "./f0fits/2mp_intens_"+str(b)+"_"+str(tBin)+".pdf"
			argandName = "./f0fits/2mp_argand_"+str(b)+"_"+str(tBin)+".pdf"
			rv.writeBinToPdf(b, stdCmd = ["", intensName, [],  argandName, []])
		return

	doF0Fits      = False
	doFits1500    = False
	doRhoFits     = False
	doRhoPfits    = False
	doRhoFfits    = False
	doF2Fits      = True
	doF2primeFits = True

	if doF2primeFits and not doF2Fits:
		raise RuntimeError("f2' fits only after f2 fits possible")
	if doF2primeFits:
		f2primeFileName = "2mp_f2primeMassesAndWidths_"+str(tBin)+".dat"
		try:
			os.remove(f2primeFileName)
		except:
			pass

	if doF0Fits:
		deg = 0
		f0fit = doF0Fit(inFileName, sectors[:1], startBin, stopBin, tBins, sectorRangeMap = {"2-+0+[pi,pi]0++PiD":(0.9,1.1)}, referenceWave = referenceWave, deg = deg, writeResultToFile = "f0MassesAndWidths_2-+0+0++_global.dat")
		with open("./f0MassesAndWidths_2mp_"+str(tBin)+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				try:
					startPars = [mF0+startValueOffset,g1+startValueOffset,g2+startValueOffset]
#					for _ in range(deg):
#						startPars.append(startValueOffset)
#						startPars.append(startValueOffset)
					x,err,c2,ndf = f0fit.fitShapeParametersForBinRange(startPars, [0],[i], zeroModeParameters = resolvedWA)
				except:
					print "Fitter exception encountered"
					startValueOffset += 0.001
					exceptCount      += 1	
					if exceptCount > 3:
						raise Exception("Too many failed attempts: "+str(exceptCount))
	
				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1])+ ' ' + str(x[2]) + ' ' + str(err[2]))
#				for i in range(deg):
#					outFile.write( ' ' + str(x[3+i]) + ' ' + str(err[3+i]) + ' ' + str(x[4+i]) + ' ' +str(err[4+i]))
				outFile.write(' ' + str(c2/ndf)+'\n')
	if doFits1500:
		with open("./2mp_f0_1500massesAndWidths_"+str(tBin)+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i + startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				if True:
					startPars = [m1500+startValueOffset,G1500+startValueOffset]
#					for _ in range(deg):
#						startPars.append(startValueOffset)
#						startPars.append(startValueOffset)
					x,err,c2,ndf = fit1500.fitShapeParametersForBinRange(startPars, [0],[i], zeroModeParameters = resolvedWA)
#				except:
#					print "Fitter exception encountered"
#					startValueOffset += 0.001
#					exceptCount      += 1
#					if exceptCount > 3:
#						raise Exception("Too many failed attempts: "+str(exceptCount))

				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
#				for i in range(deg):
#					outFile.write( ' ' + str(x[3+i]) + ' ' + str(err[3+i]) + ' ' + str(x[4+i]) + ' ' +str(err[4+i]))
				outFile.write(' ' +str(c2/ndf)+'\n')
		return


	if doRhoFits:
		with open("rhoMassesAndWidths_2mp_"+str(tBin)+rhoRangeString+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.01
				exceptCount      = 0
				while True:
					try:
#					if True:
						x,err,c2,ndf = fitBothRho.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0],[i], zeroModeParameters = resolvedWA)
						break
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							print "Too many failed attempts in bin "+str(i)+": "+str(exceptCount)
#							raise Exception
							x, err = [0.,0.],[0.,0.]
							break

				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				outFile.write(' ' + str(c2/ndf)+'\n')
	if doRhoPfits:
		with open("rhoMassesAndWidths_2mpP_"+str(tBin)+rhoRangeString+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.01
				exceptCount      = 0
				while True:
					try:
#					if True:
						x,err,c2,ndf = fitRhoP.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0],[i], zeroModeParameters = resolvedWA)
						break
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							print "Too many failed attempts in bin "+str(i)+": "+str(exceptCount)
#							raise Exception
							x, err = [0.,0.],[0.,0.]
							break

				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				outFile.write(' ' + str(c2/ndf)+'\n')

	if doRhoFfits:
		with open("rhoMassesAndWidths_2mpF_"+str(tBin)+rhoRangeString+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.01
				exceptCount      = 0
				while True:
					try:
#					if True:
						x,err,c2,ndf = fitRhoF.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0],[i], zeroModeParameters = resolvedWA)
						break
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							print "Too many failed attempts in bin "+str(i)+": "+str(exceptCount)
#							raise Exception
							x, err = [0.,0.],[0.,0.]
							break

				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				outFile.write(' ' + str(c2/ndf)+'\n')

	if doF2Fits:
		with open("f2MassesAndWidths_2mp_"+str(tBin)+f2RangeString+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.01
				exceptCount      = 0
				while True:
#					try:
					if True:
						fitF2.initMinuitFunction(fitF2.getParameters())
						x,err,c2,ndf = fitF2.fitShapeParametersForBinRange([mF2+startValueOffset,GF2+startValueOffset], [0],[i], zeroModeParameters = resolvedWA)
						break
#					except:
#						print "Fitter exception encountered"
#						startValueOffset += 0.001
#						exceptCount      += 1	
#						if exceptCount > 3:
#							print "Too many failed attempts in bin "+str(i)+": "+str(exceptCount)
##							raise Exception
#							x, err = [0.,0.],[0.,0.]
#							break

				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				outFile.write(' ' + str(c2/ndf)+'\n')

				fitF2.calculateNonShapeParametersForZeroModeParameters(resolvedWA)
				f2RV = fitF2.produceResultViewer(resolvedWA,"2-+0+[pi,pi]2++PiS", noRun = True, plotTheory = True)
				f2RV.plotData = True
				plotNameBase = "./f2FitPlots/2mp0p2ppPiS_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
				f2RV.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])


				if doF2primeFits:
					fitF2prime(i,f2primeFileName, inFileName, startBin, stopBin, tBins, [x[0], x[1]], resolvedWA, referenceWave = referenceWave,
					plotNameBase = "./f2primeFitPlots/"+str(binIndex)+"_"+str(tBin)+"_<mode>.pdf")
#					)


	if doRhoFits or doRhoPfits or doRhoFfits or doF2Fits or doF0Fits:
		return
##### Writing starts here

	fileNames = {}

	for stu in allMethods:
		print "Writing for '" + stu + "'"
		for s, sect in enumerate(allMethods[stu].sectors):
			if stu == "pipiS":
				rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, plotTheory = True)
				rv.run()
			rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, noRun = True)
			for bin in range(startBin, stopBin):
				fileName = "./collectedMethods/"+stu+"_"+sect+"_2mpData_"+str(bin)+studyAdder
				if not (sect,bin) in fileNames:
					fileNames[sect,bin] = []
				fileNames[sect,bin].append(fileName)
				rv.writeAmplFiles(bin, fileName = fileName)

	totalHists = fixedShapes.getTotalHists(resolvedWA)
	with root_open("./totals_2mp"+studyAdder+".root", "UPDATE") as out:
		for t in totalHists:
			for m in t:
				m.Write()
	folder = "./comparisonResultsData"+studyAdder+"/"

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
		allMethods['fixedShapes'].removeZeroModeFromComa()
		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = allMethods['fixedShapes'].produceResultViewer(resolvedWA,s, noRun = True, plotTheory = True)
		rv.writeBinToPdf(startBin, stdCmd = [ folder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
		for b in range(startBin, stopBin):
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			rv.writeBinToPdf(b, stdCmd = ["", folder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  folder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
	print studyList


if __name__ == "__main__":
	main()
