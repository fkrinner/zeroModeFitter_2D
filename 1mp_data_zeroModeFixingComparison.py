import matplotlib
matplotlib.use("Agg")

from   analysisClass              import amplitudeAnalysis
import parameterTrackingParameterizations as ptc
import parameterizationClasses            as pc
from   fixedparameterizationPaths import getFileNameForSector
from   modes                      import PHASE, AMPL, SMOOTH, NONE
from   utils                      import sumUp, weightedSum, cloneZeros, checkLaTeX
import sys
import pyRootPwa
import numpy as np
import numpy.linalg as la
from   globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0
from   rootfabi          import root_open
import LaTeX_strings

import consistencyUtils as cu

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter
import scipy.optimize

from random import random,randint
from math   import sin,cos

from LaTeX_strings import unCorrected_string, weightedAVG_string

def doFunctionFit(inFileName, funcs, startBin, stopBin, tBins, sectorRangeMap, referenceWave = "", writeResultToFile = None, acv = None, zeroModeParameters = None):
	sector = "1-+1+[pi,pi]1--PiP"
	fitter = amplitudeAnalysis(inFileName, [sector], {sector:funcs}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitter.loadData(loadIntegrals = True, referenceWave = referenceWave)
	fitter.finishModelSetup()
	fitter.mode = AMPL

	parameters = []
	for f in funcs:
		for p in f.returnParameters():
			parameters.append(p)

	fitter.initMinuitFunction(parameters)
	fitter.removeZeroModeFromComa()
	if acv is not None:
		fitter.addComaValueForZeroMode(acv)

	if zeroModeParameters is None:
		fitter.fitShapeParameters()
		fitter.calculateNonShapeParameters()
		chi2 = fitter.chi2
		ndf  = fitter.getNDFforMode()
	else:
		x,err,chi2,ndf = fitter.fitShapeParametersForBinRange([p.value for p in parameters], [0],range(stopBin-startBin), zeroModeParameters = zeroModeParameters)
		fitter.calculateNonShapeParametersForZeroModeParameters(zeroModeParameters)
		for i,p in enumerate(parameters):
			print p
	if writeResultToFile is not None:
		with open(writeResultToFile, 'w') as outFile:
			outFile.write("- - - - parameters - - - -\n")
			for p in parameters:
				outFile.write(str(p)+'\n')
			outFile.write(" - - - - - fit - - - - -\nchi2/NDF: "+str(chi2)+"/"+str(ndf)+"="+str(chi2/ndf))
	return fitter

def doPietarinen(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", polDeg = 1, polDeg3Pi = 0, writeResultToFile = None, acv = None, zeroModeParameters = None):
	rhoResRe  = ptc.parameter( 1.       ,  "rhoResRe"  )
	rhoResIm  = ptc.parameter( 0.       ,  "rhoResIm"  )
	rhoPoleRe = ptc.parameter( mRho**2  ,  "rhoPoleRe" )
	rhoPoleIm = ptc.parameter(-mRho*Grho,  "rhoPoleIm" )

	thresholds = [(4*mPi**2,polDeg)]
	alphaPiPi  = ptc.parameter( 0.      ,  "alphaPiPi" )
	polyCoeffs = [alphaPiPi]
	for i in range(polDeg):
		for j in range(polDeg3Pi+1):
			c = ptc.parameter( 0.,    "C"+str(i+1)+"_"+str(j)+"_PiPi" )
			polyCoeffs.append(c)
	parameters  = [rhoResRe,rhoResIm,rhoPoleRe,rhoPoleIm]
	if polDeg == 0:
		thresholds = []
	else:
		parameters += polyCoeffs

	pietarinen    = ptc.pietarinenExpansion(parameters, 1, threshsolds = thresholds, polDeg3Pi = polDeg3Pi)
	fitPietarinen = amplitudeAnalysis(inFileName, sectors, {"1-+1+[pi,pi]1--PiP":[pietarinen]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitPietarinen.loadData(loadIntegrals = True, referenceWave = referenceWave)
	fitPietarinen.finishModelSetup()
	fitPietarinen.mode = AMPL

	fitPietarinen.removeZeroModeFromComa()
	if acv is not None:
		fitPietarinen.addComaValueForZeroMode(acv)

	if zeroModeParameters is None:
		fitPietarinen.fitShapeParameters()
		fitPietarinen.calculateNonShapeParameters()
		chi2 = fitPietarinen.chi2
		ndf  = fitPietarinen.getNDFforMode()
	else:
		x,err,chi2,ndf = fitPietarinen.fitShapeParametersForBinRange([p.value for p in parameters], [0],range(stopBin-startBin), zeroModeParameters = zeroModeParameters)
		fitPietarinen.calculateNonShapeParametersForZeroModeParameters(zeroModeParameters)

	print "mRho:", rhoPoleRe.value**.5
	print "Grho:",-rhoPoleIm.value/rhoPoleRe.value**.5
	print "chi2/NDF",chi2/ndf

	if writeResultToFile is not None:
		with open(writeResultToFile, 'w') as outFile:
			outFile.write("mRho: "+str(rhoPoleRe.value**.5)+'\n')
			outFile.write("Grho: "+str(-rhoPoleIm.value/rhoPoleRe.value**.5)+'\n')
			outFile.write("- - - - parameters - - - -\n")
			for p in parameters:
				outFile.write(str(p)+'\n')
			outFile.write(" - - - - - fit - - - - -\nchi2/NDF: "+str(chi2)+"/"+str(ndf)+"="+str(chi2/ndf))
	pietarinen.writeToFile("./pietarinenFits/pietarinenAmplitude_1.5_"+str(tBins[0])+".dat", [0.139*2 + 0.001*i for i in range(2000)], [1.5])
	return fitPietarinen

def doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None, acv = None, polDeg = None):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho      = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	model    = [rho]
	if polDeg is not None:
		params = []
		sVal = 1.
		for i in range(polDeg):
			params.append(ptc.parameter( sVal, "c"+str(i+1)+"_re" ))
			sVal = 0. # Make the start value to constant 1
			params.append(ptc.parameter( 0., "c"+str(i+1)+"_im" ))
		polynome = ptc.complexPolynomial(polDeg,params)
		model.append(polynome)

	fitRho = amplitudeAnalysis(inFileName, sectors, {"1-+1+[pi,pi]1--PiP":model}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(loadIntegrals = True, referenceWave = referenceWave)
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

alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", acv = None):
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
	fixedShapes.removeZeroModeFromComa()
	if acv is not None:
		fixedShapes.addComaValueForZeroMode(acv)
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

	inFileName        = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_exotic.root"
	sectors           = ["1-+1+[pi,pi]1--PiP"]
	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	tBins             = [tBin]

	startBin          = 13
	stopBin           = 50
	if len(sys.argv) > 3:
		startBin = int(sys.argv[2])
		stopBin  = int(sys.argv[3])

	if len(sys.argv) == 4:
		startBin = int(sys.argv[2])
		stopBin  = int(sys.argv[3])

	acv               = 130. # artificial coma value
#	acv               = 1.3e-5 # artificial coma value

	seedint = randint(0,10000)

#	startBin          = 25
#	stopBin           = 30

	methodBinRanges   = {
	                   "fitF2"        : (22, 50),
	                   "fitRhoF"      : (22, 50),
	                   "fixedShapeF2" : (22, 50)}
#	methodBinRanges   = {} # Override here 

	rhoRange = None
	for a in sys.argv:
		if a.startswith("range"):
			rhoRange = float(a[5:])
			sectorRangeMap = {"1++0+[pi,pi]1--PiS":(0.,rhoRange)}

	sectorRangeMap    = {}
#	sectorRangeMap    = {"1-+1+[pi,pi]1--PiP":(0.,1.2)}
	fixSectorRangeMap = {"1-+1+[pi,pi]1--PiP":(0.,1.12)}	

	allMethods        = {}
	methodStrings     = {}
	shortlabels       = { "fixedShapeF0"      : r"$\text{fix}_{f_0}^{~}$",
	                      "fixedShapeRhoP"    : r"$\text{fix}_\rho^{P}$",
	                      "fixedShapeRhoF"    : r"$\text{fix}_\rho^{F}$",
	                      "fixedShapeBothRho" : r"$\text{fix}_\rho^{2}$",
	                      "fixedShapeF2"      : r"$\text{fix}_{f_2}$",
	                      "fixedShapes"       : r"$\text{fix}_\text{all}^{~}$",
	                      "pipiS"             : r"$\phi_{[\pi\pi]_S}^{~}$",
	                      "fitRho"            : r"$\text{fit}_\rho$",
	                      "fitRhoP"           : r"$\text{fit}_\rho^{P}$",
	                      "fitRhoF"           : r"$\text{fit}_\rho^{F}$",
	                      "fitBothRho"        : r"$\text{fit}_\rho^{2}$",
	                      "fitF2"             : r"$\text{fit}_{f_2}$",
	                      "smooth"            : r"smooth"}

#       # - - - - --- Start here with the model builtind --- - - - - #       #
	rhoRe0     =  0.579053
	rhoIm0     =  0.109177

	rhoMass    = ptc.parameter(mRho, "rho_mass")
	rhoWidth   = ptc.parameter(Grho, "rho_width")
	rho        = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)

	poleReal   = ptc.parameter(rhoRe0, "rhoRe")
	poleImag   = ptc.parameter(rhoIm0, "rhoIm")

	nPol = 1
	for a in sys.argv:
		if a.startswith("nPol"):
			nPol = int(a[4:])
	polyDeg_po = 8
	for a in sys.argv:
		if a.startswith("pdpo"):
			polyDeg_po = int(a[4:])


	params     = [poleReal,poleImag]
	for i in range(1, nPol):
		rePol = ptc.parameter(1.96, "polRe_"+str(i))
		imPol = ptc.parameter(0.28, "polIm_"+str(i))
		params.append(rePol)
		params.append(imPol)

	for d in range(polyDeg_po):
		params.append(ptc.parameter(2.*random()-1., "c_"+str(d)))
	Kmatrix    = ptc.simpleOneChannelKmatrix(params, nPol, polyDeg_po, 4*mPi**2)
	useCM      = False
	for a in sys.argv:
		if a == "CM":
			useCM = True
		if a == "rho":
			useCM = False

	Kmatrix.use_CM = useCM

	pPolyDeg3  = 7
	pPolyDeg2  = 4
	for a in sys.argv:
		if a.startswith("m3pol"):
			pPolyDeg3 = int(a[5:])
		if a.startswith("m2pol"):
			pPolyDeg2 = int(a[5:])


	params     = []
	for d in range(pPolyDeg2):
		for e in range(pPolyDeg3+1):
			params.append(ptc.parameter(2.*random()-1., "c_"+str(d+1)+"_"+str(e)))
	pPoly      = ptc.twoDimensionalRealPolynomial(pPolyDeg2, pPolyDeg3, params, baseExponent = 2) # baseExponent = 2: polynomial in s
	func       = ptc.multiply([Kmatrix, pPoly])	

	model      = [func]
#	model      = [rho]

#       # - - - - --- Stop the model building here --- - - - - #       #
	zeroModeParameters = None
	fixZeroMode        = True
	if fixZeroMode:
		fixedShapes        = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, acv = acv, sectorRangeMap = fixSectorRangeMap)
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
	resultFile  =  "./KmatrixResults/1mp_Kmatrix_"+rrs+nps+ps+"_kPol"+str(polyDeg_po-1)+"_pPol"+str(pPolyDeg3)+"-"+str(pPolyDeg2)+"_t"+str(tBin)+"_m"+str(startBin)+'-'+str(stopBin)+'_'+str(seedint)+".dat"
	fitter      = doFunctionFit(inFileName, model, startBin, stopBin, tBins, sectorRangeMap,    referenceWave = referenceWave, acv = acv, zeroModeParameters = zeroModeParameters, writeResultToFile = resultFile)
	rhoFitter   = doFunctionFit(inFileName, [rho], startBin, stopBin, tBins, fixSectorRangeMap, referenceWave = referenceWave, acv = acv, zeroModeParameters = zeroModeParameters, writeResultToFile = None)
	BWfitString = "BW_fit_result: " + str(rhoMass.value) + " " +str(rhoMass.error) + " " + str(rhoWidth.value) + " " + str(rhoWidth.error)

#       # - - - - --- Start the evaluations here --- - - - - #       #
	nBinsPlot = 1000
	def fPlot(v):
		return v.imag
#	hist = pyRootPwa.ROOT.TH2D("hhh","hhh", nBinsPlot, -.25, 6.25,nBinsPlot, -1., 1.)
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
	res = scipy.optimize.minimize(Kmatrix.absInverse,[rhoRe0,rhoIm0])
	print res.x,"pole position"
	mfv = res.fun
	resSting = str(res.fun)+ " function value should be zero"
	print resSting
	BWstring = "BW par: "+str(abs(res.x[0])**.5)+" "+str(abs(res.x[1])/abs(res.x[0])**.5)+" (all absolute values)"
	print BWstring


	with open(resultFile,'a') as outFile:
		outFile.write('\n'+BWstring+" "+resSting)

	print "= = = = = = = = = Starting BW error ersimation = = = = = = = = = "
	nPoints     = 1000
	poleMean    = res.x
	poleSamples = []
	i           = 0
	failCount   = 0
	while i < nPoints:
		pts = np.random.multivariate_normal(fitter.fitParameters, fitter.MINUITcoma)
		fitter.MINUIT_function(pts) # Call the function once to set parameters inside
#		fitter.model[0].setParametersAndErrors(pts, fitter.MINUITerrs)
		res = scipy.optimize.minimize(Kmatrix.absInverse,poleMean)
		if abs(res.fun) > 100*mfv:
			print "No more pole found (mfv = "+str(mfv)+") : fval = "+str(res.fun)
			failCount += 1
			if failCount > nPoints:
				print "Failed to find poles too often.... abort"
				return
			continue
#			raise ValueError("No more pole found: fval = "+str(res.fun))

		poleSamples.append(res.x)
		i+=1
		print i,"Marker to find the PRINT 57473M3N7"
	meanPole = [0.,0.]
	for p in poleSamples:
		meanPole[0] += p[0]
		meanPole[1] += p[1]
	meanPole[0] /= len(poleSamples)
	meanPole[1] /= len(poleSamples)
	poleComa = [[0.,0.],[0.,0.]]
	for p in poleSamples:
		poleComa[0][0] += (p[0]-meanPole[0])**2
		poleComa[0][1] += (p[0]-meanPole[0])*(p[1]-meanPole[1])
		poleComa[1][0] += (p[1]-meanPole[1])*(p[0]-meanPole[0])
		poleComa[1][1] += (p[1]-meanPole[1])**2
	poleComa[0][0] /= len(poleSamples)-1
	poleComa[0][1] /= len(poleSamples)-1
	poleComa[1][0] /= len(poleSamples)-1
	poleComa[1][1] /= len(poleSamples)-1
	comaString      = str(poleComa)
	print " - - - - - - le compaire pramaitre  - - - - - - "
	print meanPole, poleMean
	print " - - - - - - le compaire pramaitre  - - - - - - "
	print poleComa
	print "= = = = = = = = = Finished BW error ersimation = = = = = = = = = "
	mRhoP = 1.4
	GrhoP = .2
	res = scipy.optimize.minimize(Kmatrix.absInverse,[mRhoP**2,mRhoP*GrhoP])
	print res.x,"pole position"
	resSting = str(res.fun)+ " function value should be zero"
	print resSting
	BWstring = "BW' par: "+str(abs(res.x[0])**.5)+" "+str(abs(res.x[1])/abs(res.x[0])**.5)+" (all absolute values)"
	print BWstring
	with open(resultFile,'a') as outFile:
		outFile.write('\n'+BWstring+" "+resSting+"\ncoma "+comaString)
		outFile.write('\n'+BWfitString)
#       # - - - - --- Start the evaluations here --- - - - - #       #
	doPlots = False
	if doPlots:
		RV = fitter.produceResultViewer(zeroModeParameters,"1-+1+[pi,pi]1--PiP", noRun = True, plotTheory = True)
		RV.plotData = True
		for b in range(startBin, stopBin):
			plotNameBase = "./Kmatrix_plots/1mp1p1mmPiP_<mode>_"+str(b)+"_"+str(tBin)+".pdf"
			RV.writeBinToPdf(b, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
#	raw_input("press <enter> to exit")
	return
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
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### #
#################################################################################
#################################################################################
# The stuff down here does not matter at the moment #
	
	polDeg      = 2
	outFileName = "./pietarinenFits/polyDeg"+str(polDeg)+"_t"+str(tBin)+".txt"

	piet = doPietarinen(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, polDeg = polDeg, polDeg3Pi = 2, writeResultToFile = outFileName, sectorRangeMap = sectorRangeMap, acv = acv, zeroModeParameters = zeroModeParameters)
	if not fixZeroMode:
		zeroModeParameters = piet.getZeroModeParametersForMode()
	RV = piet.produceResultViewer(zeroModeParameters,"1-+1+[pi,pi]1--PiP", noRun = True, plotTheory = True)
	RV.plotData = True
	for b in range(startBin, stopBin):
		plotNameBase = "./pietarinenFits/1mp1p1mmPiP_<mode>_"+str(b)+"_"+str(tBin)+".pdf"
		RV.writeBinToPdf(b, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
#	doFitRho(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "rhoMassesAndWidths_1-+1+1--_global.dat", sectorRangeMap = sectorRangeMap, acv = acv, polDeg = polDeg)
	return

	print "Starting with fixed shapes"
	fixedShapes = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, acv = acv, sectorRangeMap = sectorRangeMap)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"

#	print "Starting with fitting rho"
#	fitRho = doFitRho(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "rhoMassesAndWidths_1-+1+1--_global.dat", sectorRangeMap = sectorRangeMap, acv = acv)
#	allMethods["fitRho"] = fitRho
#	print "Finished with fitting rho"

#	if stopBin - startBin > 1:
#		print "Starting with smooth"
#		smooth = doSmooth(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
#		allMethods["smooth"] = smooth
#		print "Finished with smooth"

	if "fixedShapeRhoF" in allMethods:
		allMethods["fixedShapeRhoF"].setZeroModeSignature(fullSig,1)
	if "fitRhoF" in allMethods:
		allMethods["fitRhoF"].setZeroModeSignature(fullSig,1)

#	diffsFull, resolvedWA, nonResolvedWA, comps, resDiffs, nonResDiffs, resolvedDiffsFull,noCorrDiffs = cu.doAllComparisons(allMethods, startBin, methodBinRanges, noBelowZero = True)

#	print resolvedDiffsFull
#	from math import isnan
#	for pair in resolvedDiffsFull:
#		with  modernplotting.toolkit.PdfWriter('./resolvedDiffPlots/1mp_'+pair[0]+"_"+pair[1]+"_"+str(tBin)+".pdf") as pdfOutput:
#			plot  = style.getPlot1D()
#			line  = [0.000000001]*len(resolvedDiffsFull[pair][0])
#			line2 = [0.000000001]*len(resolvedDiffsFull[pair][0])
#			one   = [1.]*len(resolvedDiffsFull[pair][0])
#			xAxis = [ .5 + 0.04*(startBin + i) for i in range(len(resolvedDiffsFull[pair][0]))]
#			for i,v in enumerate(resolvedDiffsFull[pair][0]):
#				if not isnan(v) and not v <= 0.:
#					line[i] = v
#				else:
#					line[i] = 0.000000001
#
#			if not pair[1] == "WAres" and not pair[1] == "WAnon":
#				for i,v in enumerate(resolvedDiffsFull[pair[1],pair[0]][0]):
#					if not isnan(v) and not v <= 0.:
#						line2[i] = v
#					else:
#						line2[i] = 0.000000001
#
#			plot.setYlog()
#			plot.plot(xAxis, line)
#			plot.plot(xAxis, one)
#			plot.plot(xAxis, line2)
#			plot.setYlim(0.00001, 10000)
#			plot.finishAndSaveAndClose(pdfOutput)

#	studyList = []
#	for m in allMethods:
#		studyList.append(m)
#	studyList.sort()

	style.titleRight = r"$1^{-+}1^+$"
	style.titleLeft  = LaTeX_strings.tBins[tBin]

#	with  modernplotting.toolkit.PdfWriter("compositions_1mp_"+str(tBin)+".pdf") as pdfOutput:
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
#		plot.finishAndSaveAndClose(pdfOutput)

#	hist = pyRootPwa.ROOT.TH2D("hist","hist", len(studyList)+2, 0, len(studyList)+2, len(studyList), 0, len(studyList))
#
#	for i,m in enumerate(studyList):
#		for j,n in enumerate(studyList):
#			hist.SetBinContent(i+1, j+1, diffsFull[n,m])
#	for i,m in enumerate(studyList):
#		hist.SetBinContent(len(studyList)+1, i+1, noCorrDiffs[m])
#		hist.SetBinContent(len(studyList)+2, i+1, resDiffs[m])
#
#	axolotl = []
#	for i,study in enumerate(studyList):
#		axolotl.append(shortlabels[study])
#		axolotl.append(alphabet[i])

#	with modernplotting.toolkit.PdfWriter("studies_1mp_"+str(tBin)+".pdf") as pdfOutput:
#		plot = style.getPlot2D()
#		plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(len(studyList)+2)])
#		plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(len(studyList)  )])
#		studyPlotter.makeValuePlot(plot, hist)
#		
#		plot.axes.set_yticklabels(axolotl)
#		axolotl.append(unCorrected_string)
#		axolotl.append(weightedAVG_string)
#		plot.axes.set_xticklabels(axolotl, rotation = 90)
#		plot.setZlim((0.,1.))
#
#		plot.finishAndSaveAndClose(pdfOutput)

#	with open("studies_1mp_"+str(tBin)+".txt",'w') as out:
#		for axl in axolotl:
#			out.write(axl + ' ')
#		out.write("\n")
#		for i in range(hist.GetNbinsX()):
#			for j in range(hist.GetNbinsY()):
#				out.write(str(hist.GetBinContent(i+1, j+1)) + ' ')
#			out.write('\n')

	doRhoFits = False
	writeCpls = False
	if writeCpls:
		outFileCpl = open("1mp_rho_cpls_"+str(tBin)+".dat",'w') 

	doActuallyNotFit = True
	if doRhoFits:
		with open("rhoMassesAndWidths_exotic_"+str(tBin)+".dat",'w') as out:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				out.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				if doActuallyNotFit:
					print "The fit has been turned off, due to a workaround... :("
				else:
					exceptCount      = 0
					try:
						x,err,c2,ndf = fitRho.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0],[i], zeroModeParameters = resolvedWA)
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							raise Exception("Too many failed attempts: "+str(exceptCount))
	
					out.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
					out.write(' '+str(c2/ndf)+'\n')

				if writeCpls:
					fitRho.calculateNonShapeParametersForZeroModeParameters(resolvedWA)
					cpl, hess = fitRho.getCouplingParameters()
					hessInv = la.inv(hess[0][i])
					if not len(cpl[0][1]) == 2:
						raise IndexError("Parameter count not 2, change implementation")
					integral = fitRho.model[0][i].getIntegralForFunction(0, fitRho.model[0][i].funcs[0][0])
					outFileCpl.write(str(0.52 + binIndex*.04) + ' ' + str(cpl[0][i][0]**2 + cpl[0][i][0]**2) + ' ' + str(integral) + ' ')
					outFileCpl.write(str(cpl[0][i][0])     + ' ' + str(cpl[0][i][1])     + ' ')
					outFileCpl.write(str(hessInv[0][0]/2) + ' ' + str(hessInv[1][1]/2) + ' ' + str(hessInv[0][1]/2))
					outFileCpl.write("\n")

				makeRhoFitPlots = False
				if makeRhoFitPlots:
					fitRho.calculateNonShapeParametersForZeroModeParameters(resolvedWA)
					rhoRV = fitRho.produceResultViewer(resolvedWA,"1-+1+[pi,pi]1--PiP", noRun = True, plotTheory = True)
					rhoRV.plotData = True
					plotNameBase = "./rhoFitPlots/1mp1p1mmPiP_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
					rhoRV.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
		if writeCpls:
			outFileCpl.close()

		return
##### Writing starts here

#	fileNames = {}

#	fixedShapes.removeZeroModeFromComa()
##	binRange = {"1-+1+[pi,pi]1--PiP": (10,22)}
#	totalHists = fixedShapes.getTotalHists(resolvedWA, binRange = binRange)
#	with root_open("./totals_1mp_rhoRange.root", "UPDATE") as out:
#		for t in totalHists:
#			for m in t:
#				m.Write()
#	return

#	for stu in allMethods:
#		print "Writing for '" + stu + "'"
#		for s, sect in enumerate(allMethods[stu].sectors):
#			if stu == "pipiS":
#				rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, plotTheory = True)
#				rv.run()
#			rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, noRun = True)
#			for bin in range(startBin, stopBin):
#				fileName = "./collectedMethods/"+stu+"_"+sect+"_1mpData_"+str(bin)
#				if not (sect,bin) in fileNames:
#					fileNames[sect,bin] = []
#				fileNames[sect,bin].append(fileName)
#				rv.writeAmplFiles(bin, fileName = fileName)
	scalle = False
	if scalle:
		folder = "./comparisonResultsData_1mp_scale/"
	else:
		folder = "./comparisonResultsData_1mp/"

	fixedShapesFull = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, acv = acv)

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
#		allMethods['fixedShapes'].removeZeroModeFromComa()
#		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = fixedShapesFull.produceResultViewer(allMethods['fixedShapes'].getZeroModeParametersForMode(),s, noRun = True, plotTheory = True)
#		rv.plotData = False
		rv.writeBinToPdf(startBin, stdCmd = [ folder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
#		rv.labelPoints     = [0,10,15,20]
#		rv.makeLegend      = True
		if scalle:
			rv.scaleTo         = "maxCorrData"
#			rv.yAxisShift      = 300.
#			rv.tStringYpos     = 0.8
#			rv.topMarginIntens = 1.4
#			fakkkk             = 1.
		else:
			rv.plotData        = False
#			fakkkk             = .7	
#			rv.tStringYpos     = 0.865
#			rv.topMarginIntens = 1.3
#			rv.yAxisShift      = 100.
#		rv.shiftMap        = {0:(fakkkk*50.,fakkkk*-280.),10:(fakkkk*-420.,fakkkk*-50.), 15:(fakkkk*-420.,fakkkk*-30.), 20:(fakkkk*-50.,fakkkk*70.)}
		for b in range(startBin, stopBin):
#			if not bin == 27:
#				continue
#			intensNames = [name+".intens" for name in fileNames[sect,b]]
#			argandNames = [name+".argand" for name in fileNames[sect,b]]
#			intensNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/evalDima/dima.intens"]
#			argandNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/evalDima/dima.argand"]
#			if not scalle:
			if True:
				intensNames = []
				argandNames = []
#			rv.plotData = True
			rv.writeBinToPdf(b, stdCmd = ["", folder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  folder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
			rv.wiriteReImToPdf(b, folder +  sect + "_data_<ri>_"+str(b)+"_"+str(tBin)+".pdf")


if __name__ == "__main__":
	main()
