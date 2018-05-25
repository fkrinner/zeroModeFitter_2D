import matplotlib
matplotlib.use("Agg")

from analysisClass import amplitudeAnalysis
import parameterTrackingParameterizations as ptc
import parameterizationClasses as pc
from fixedparameterizationPaths import getFileNameForSector
from modes import PHASE, AMPL, SMOOTH, NONE
from utils import sumUp, weightedSum, cloneZeros, checkLaTeX
import sys
import pyRootPwa
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
	sector  =  "2-+0+[pi,pi]2++PiS"
	sectors = ["2-+0+[pi,pi]0++PiD","2-+0+[pi,pi]1--PiP","2-+0+[pi,pi]1--PiF","2-+0+[pi,pi]2++PiS"]

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

def main():
	checkLaTeX()
	style = modernplotting.mpplot.PlotterStyle()
#	style.p2dColorMap = 'ocean_r'
#	style.p2dColorMap = 'YlOrRd'
	style.p2dColorMap = 'Reds'
#	referenceWave     = ""

	inFileName = fileNameMap["std11"]

	for a in sys.argv:
		if a.startswith("mBins"):
			startBin = int(a.split("-")[1])
			stopBin  = int(a.split("-")[2])
		if a.startswith("tBins"):
			tBins = range(int(a.split("-")[1]),int(a.split("-")[2]))


	nPol       = 3
	for a in sys.argv:
		if a.startswith("nPol"):
			nPol = int(a[4:])

	f2Re0     = mF2**2
	f2Im0     = mF2*GF2

	mPrime = 1.55
	Gprime =  .2

	mPPrime = 1.91
	Gpprime = .3


	f2Mass     = ptc.parameter(mF2, "f2_mass")
	f2Width    = ptc.parameter(GF2, "f2_width")
	f2         = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, False)

	poleReal   = ptc.parameter(f2Re0, "f2Re")
	poleImag   = ptc.parameter(f2Im0, "f2Im")

	pPoleReal  = ptc.parameter(mPrime**2,      "f2pRe")
	pPoleImag  = ptc.parameter(mPrime*Gprime, "f2pIm")

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

	for d in range(polyDeg_po):
		params.append(ptc.parameter(2*random()-1., "c_"+str(d)))

	Kmatrix    = ptc.simpleOneChannelKmatrix(params, nPol, polyDeg_po, 4*mPi**2)
	useCM      = False
	for a in sys.argv:
		if a == "CM":
			useCM = True
		if a == "rho":
			useCM = False
	Kmatrix.use_CM = useCM

	for a in sys.argv:
		if a.startswith("polyDegP"):
			polyDegP = int(a[8:])

	mMin = .5 + .04*startBin
	mMax = .5 + .04*stopBin

	tFuncs = []
	tpPoly = []
	nKpar  = len(params)

	for tBin in tBins:
		tParams = []
		for mBin in range(startBin,stopBin):
			for d in range(polyDegP):
				tParams.append(ptc.parameter(2*random()-1., "c"+str(d+1)+"_t"+str(tBin)+"_m"+str(mBin)))
		nPpar   = len(tParams)
		params += tParams
		pPoly   = ptc.binnedPolynomial(tParams, mMin = mMin, mMax = mMax, nBins = stopBin-startBin, degree = polyDegP, baseExponent = 2, real = True)
		tpPoly.append(pPoly)
		func    = ptc.multiply([Kmatrix, pPoly])
		tFuncs.append(func)

	sector  =  "2-+0+[pi,pi]2++PiS"
	sectors = ["2-+0+[pi,pi]0++PiD","2-+0+[pi,pi]1--PiP","2-+0+[pi,pi]1--PiF","2-+0+[pi,pi]2++PiS"]

	sectorRangeMap = {}
	zmPars         = []
	fitters        = []
	for t,tBin in enumerate(tBins):
		fitter = amplitudeAnalysis(inFileName, sectors, {sector:[tFuncs[t]]}, startBin, stopBin, [tBin], sectorRangeMap = sectorRangeMap)
		fitter.loadData(loadIntegrals = True, referenceWave = referenceWave)
		fitter.finishModelSetup()
		fitter.mode = AMPL

		fixSectors = ["2-+0+[pi,pi]1--PiP", "2-+0+[pi,pi]1--PiF", "2-+0+[pi,pi]2++PiS"]
		fixRangeMap = {
			"2-+0+[pi,pi]1--PiP" : (0.,1.12),
			"2-+0+[pi,pi]1--PiF" : (0.,1.12),
			"2-+0+[pi,pi]2++PiS" : (0.,1.27) # Use only left side of f_2(1270)
		}

		fixedShapes        = doFixedShapes(inFileName, fixSectors, startBin, stopBin, [tBin], referenceWave = referenceWave, sectorRangeMap = fixRangeMap)
		fitter.setZeroModeParameters(fixedShapes.getZeroModeParametersForMode())
		fitter.model.setBinsToEvalueate([0], range(stopBin-startBin)) # [0] since it is only one tBin in this fitter (one fitter par t' bin)
		fitters.append(fitter)

	def evalFunction(params):
		Kmatrix.setParameters(params[:nKpar])
		for t in range(len(tBins)):
			tpPoly[t].setParameters(params[nKpar*t*nPpar:nKpar*(t+1)*nPpar])
		chi2 = 0.	
		for t in range(len(tBins)):
			chi2 += fitters[i].fixedZMPchi2(None)

	print evalFunction([p.value for p in params]), "gande gange gand gond gale"

	return

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
	if nPol == 1:
		nps = ""
	else:
		nps = "nPol"+str(nPol)+"_"

	resultFile  =  "./KmatrixRestrictedZMfixing_allFourSectorsActive/2mpF2_Kmatrix_"+nps+ps+"_kPol"+str(polyDeg_po-1)+"_pPol"+str(pPolyDeg3)+"-"+str(pPolyDeg2)+"_t"+str(tBin)+"_m"+str(startBin)+'-'+str(stopBin)+'_'+str(seedint)+".dat"
	fitter      = doFunctionFit(inFileName, model, startBin, stopBin, tBins, sectorRangeMap, referenceWave = referenceWave, acv = acv, zeroModeParameters = zeroModeParameters, writeResultToFile = resultFile, ifn = ifn)

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
	res = scipy.optimize.minimize(Kmatrix.absInverse,[f2Re0,f2Im0])
	print res.x,"pole position"
	mfv = res.fun
	resSting = str(res.fun)+ " function value should be zero"
	print resSting
	BWstring = "BW par: "+str(abs(res.x[0])**.5)+" "+str(abs(res.x[1])/abs(res.x[0])**.5)+" (all absolute values)"
	print BWstring

	if ifn is None:
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
	#		print i,"Marker to find the PRINT 57473M3N7"
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
		mF2P = 1.9
		GF2P =  .277
		res  = scipy.optimize.minimize(Kmatrix.absInverse,[mF2P**2,mF2P*GF2P])
		print res.x,"pole position"
		resSting = str(res.fun)+ " function value should be zero"
		print resSting
		BWstring = "BW' par: "+str(abs(res.x[0])**.5)+" "+str(abs(res.x[1])/abs(res.x[0])**.5)+" (all absolute values)"
		print BWstring
		with open(resultFile,'a') as outFile:
			outFile.write('\n'+BWstring+" "+resSting+"\ncoma "+comaString)
	#		outFile.write('\n'+BWfitString)
	#       # - - - - --- Start the evaluations here --- - - - - #       #
		doPlots = False
		for a in sys.argv:
			if a == "plot":
				doPlots = True
	if ifn is not None:
		doPlots = True # Set plots by default, if an inFile is given

	if doPlots:
		RV = fitter.produceResultViewer(zeroModeParameters,"2-+0+[pi,pi]2++PiS", noRun = True, plotTheory = True)
		RV.plotData = True
		for b in range(startBin, stopBin):
			plotNameBase = "./Kmatrix_plots/2mp0p2ppPiS_<mode>_"+str(b)+"_"+str(tBin)+"_"+str(seedint)+".pdf"
			RV.writeBinToPdf(b, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
#	raw_input("press <enter> to exit")
	return

if __name__ == "__main__":
	main()
