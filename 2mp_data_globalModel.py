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
from cmath import pi

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from random import randint, random
import scipy
from iminuit import Minuit
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

def explicit_MINUIT_function():
	pass

def function():
	pass

def MIGRAD(argString, lstString, function2):
	global function
	function = function2
	command = "def eMf("+lstString+"):\n\treturn function(["+lstString+"])"
	global explicit_MINUIT_function
	exec command
	explicit_MINUIT_function = eMf

	command ="m = Minuit(explicit_MINUIT_function,"+argString+")"
	exec command
	print "Starting migrad"
	m.migrad()

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
			tMin  = int(a.split("-")[1])
			tMax  = int(a.split("-")[2])
			tBins = range(tMin, tMax)

	nPol       = 3
	for a in sys.argv:
		if a.startswith("nPol"):
			nPol = int(a[4:])

	f2Re0   = mF2**2
	f2Im0   = mF2*GF2

	mPrime  = 1.55
	Gprime  =  .2

	mPPrime = 1.91
	Gpprime =  .3

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
	if useCM:
		psString = "CM"
	else:
		psString = "rho"

	Kmatrix.use_CM = useCM

	for a in sys.argv:
		if a.startswith("polyDegP"):
			polyDegP = int(a[8:])

	mMin = .5 + .04*startBin
	mMax = .5 + .04*stopBin

	tFuncs = []
	tpPoly = []
	nKpar  = len(params)
	KparValues = [p.value for p in params]
	nKpar  = 0

	funcs = []
	for d in range(polyDegP+1):
		mnm  = ptc.monomial(2*d)
		func = ptc.multiply([Kmatrix, mnm])
		funcs.append(func)
	bfv = [	 1.591236,   0.234037,   2.419353,   0.005742,   3.702362,   0.000053,  -0.189620, 0.207497 ,-0.058921,0.517699,0.166935 ,0.181275,0.959202 , 1.528309, 0.501027, -1.146739,-0.535555,  0.417481, -0.419879,  1.561247,
		 0.120558,   1.059373,   0.393481,  -0.283450,  -0.182751,  -0.667384,   0.619867,   0.343112, 0.069652 , 0.677225, -0.116518, -0.485594,  0.693124,  0.270239,  0.139340, -2.694811, -0.062632,  1.570524, 
		-0.493674,  -2.409801,   0.394341,   0.421391,  -0.524359,   0.112583,   0.285943,  -0.166501, -0.408990, -0.010904,  0.088735,  0.117378,  0.166603, -0.949778,  0.500537, -1.352419, -0.639588, -2.139980, 
		-0.146114,   1.877676,   1.025465,  -0.120469,   1.197397,   0.792957,  -0.507836,   0.354736, -0.571950, -0.439970,  0.635951,  1.119741,  1.078536, -1.707963, -0.894598, -0.379259, -0.256336,  1.811274,  
		 0.460570,  -2.370441,  -1.519688,  -2.596059,  -1.472923,   0.329489,   0.963846,  -0.174265, 0.497754,  0.523653, -0.082330, -2.288272,  0.238584,  2.432317, -0.219072, -1.818075,  0.101130, -1.315669, 
		-0.095197,   0.185721,  -0.468701,   1.150421,  -0.706067,   1.336765,   0.114122,  -2.514466, -1.292956,  0.909066, -0.136537, 1.579291,  0.917968,  0.519039,  0.560456,  -3.141557,  -0.919017,   1.034258,   
		 0.651617,   0.414377,  -0.954656,   0.432315,  -0.250071,  -0.306383,  -0.916611,   0.241417,   0.196445,  -0.594122,  -0.610276,  -2.115296,  -1.232504,   0.210273,  -0.539883,   1.134149,   0.582253,  
		-1.330798,  -1.189863,   2.000385,  -0.286371,   0.997649,  -0.268911,  -0.358648,   0.209653,   0.307646,   0.294261,   2.350279,  -0.342739,   1.010674,  -0.547855,  -1.715226,   0.599891,   0.567254,
		-0.731444,   2.852883,  -0.398756,   2.450370,   0.296712,   0.499203,   0.284769,  -7.040768,  -1.148202,  -0.982813,   0.349752,  -6.325548,  -0.494806,   0.183253,   0.745206,   2.547205,   1.571422,
		 0.713172,  -0.253361,   0.610771,   1.308658,   1.501376,   1.023860,   2.184354,   0.927320,   0.677329,   0.589890,   0.139710,   0.550017,  -0.702369,  -0.847349,  -0.729193,   0.541576,  -0.508911, 
		 0.495425,   1.287986,  -0.364052,  -1.338609,   0.281415,   0.168475,   0.184298,  -0.904410,   0.193151,  -0.773044,   0.324601,   4.473005,   0.046860,  -0.359312,   0.368391,  -1.623123,   0.531775, 
		 0.384091,   0.988951,  -0.779931,   0.030888,  -0.139630,  -0.264398,  -9.821342,  -1.601821,  -7.745628,  -2.010218,  -3.540738,  -0.067302,   0.153977,-1.635704]
#	bestFitVals = [1.591236, 0.234037,2.419353 ,0.005742, 3.702362, 0.000053 , -0.189620 ,0.207497 ,-0.058921]
	for p,par in enumerate(params):
		par.value = bestFitVals[p]

	phases = []
	for tBin in tBins:
#		tParams = []
		for mBin in range(startBin,stopBin):
			phases.append(ptc.parameter(2*pi*random(), "phi_t"+str(tBin)+"_m"+str(mBin)))
#				tParams.append(ptc.parameter(2*random()-1., "c"+str(d+1)+"_t"+str(tBin)+"_m"+str(mBin)))
#		nPpar   = len(tParams)
#		params += tParams
#		pPoly   = ptc.binnedPolynomial(tParams, mMin = mMin, mMax = mMax, nBins = stopBin-startBin, degree = polyDegP, baseExponent = 2, real = True)
#		tpPoly.append(pPoly)
#		func    = ptc.multiply([Kmatrix, pPoly])
#		tFuncs.append(func)

#	params += phases
	params = phases

	sector  =  "2-+0+[pi,pi]2++PiS"
	sectors = ["2-+0+[pi,pi]0++PiD","2-+0+[pi,pi]1--PiP","2-+0+[pi,pi]1--PiF","2-+0+[pi,pi]2++PiS"]

	sectorRangeMap = {}
	zmPars         = []
	fitters        = []
	for t,tBin in enumerate(tBins):
		fitter = amplitudeAnalysis(inFileName, sectors, {sector:funcs}, startBin, stopBin, [tBin], sectorRangeMap = sectorRangeMap)
		fitter.loadData(loadIntegrals = True, referenceWave = referenceWave)
		fitter.finishModelSetup()
		fitter.mode = AMPL

		fixSectors = ["2-+0+[pi,pi]1--PiP", "2-+0+[pi,pi]1--PiF", "2-+0+[pi,pi]2++PiS"]
		fixRangeMap = {
			"2-+0+[pi,pi]1--PiP" : (0.,1.12), 			"2-+0+[pi,pi]1--PiF" : (0.,1.12), 			"2-+0+[pi,pi]2++PiS" : (0.,1.27) # Use only left side of f_2(1270)
		}

		fixedShapes = doFixedShapes(inFileName, fixSectors, startBin, stopBin, [tBin], referenceWave = referenceWave, sectorRangeMap = fixRangeMap)
		fitter.setZeroModeParameters(fixedShapes.getZeroModeParametersForMode())
		fitter.model.setBinsToEvalueate([0], range(stopBin-startBin)) # [0] since it is only one tBin in this fitter (one fitter par t' bin)
		fitters.append(fitter)

	nmBins = stopBin - startBin
	def evalFunction(params):
#		Kmatrix.setParameters(params[:nKpar])
		chi2 = 0.	
		for t in range(len(tBins)):
#			chi2 += fitters[t].model.fixedZMPchi2(None)
			chi2 += fitters[t].model.fixedZMPchi2_realCouplings(None, params[nKpar+t*nmBins:nKpar+(t+1)*nmBins])
		return chi2

	argString = ','.join([p.name+'='+str(p.value) for p in params])
	lstString = ','.join([p.name for p in params])

	ndf  = 0
	for fitter in fitters:
		ndfs = fitter.model.getNDF()
		for t in ndfs:
			for m in t:
				ndf += m[0] - m[1] - m[2]
	ndf -= len(params)
	ndf += polyDegP*(stopBin-startBin)*(tMax-tBin) # the couplings have jus a global phase !!!
# # # # # # 2mpF2_Kmatrix_nPol3_rho_kPol0_pPol0-3_t3_m49-50_5854.dat

	for att in range(10): # make 10 attempts at once
		print "Setting Start Pars for attempt",att+1
		for i in range(len(params)):
			if i < nKpar:
				params[i].value = KparValues[i]
			else:
				params[i].value = 2*pi*random()
		MIGRAD(argString, lstString, evalFunction)

		chi2 = explicit_MINUIT_function(*[p.value for p in params])

		seedint = randint(0,10000)
		resultFileName = "./globalKmatrixFits/2mpF2_nPol"+str(nPol)+"_"+psString+"_kPol"+str(polyDeg_po-1) + "_t"+str(tMin)+'-'+str(tMax)+"_m"+str(startBin)+"-"+str(stopBin)+"_polyDegP"+str(polyDegP)+"_"+str(seedint)+".dat"

		with open(resultFileName, 'w') as outFile:
			outFile.write("- - - - parameters - - - -\n")
			for p in params:
				outFile.write(str(p)+'\n')
			outFile.write(" - - - - - fit - - - - -\nchi2/NDF: "+str(chi2)+"/"+str(ndf)+"="+str(chi2/ndf)+'\n')
			outFile.write("ndf corrected!!!")
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
			"2-+0+[pi,pi]1--PiP" : (0.,1.12), 			"2-+0+[pi,pi]1--PiF" : (0.,1.12), 			"2-+0+[pi,pi]2++PiS" : (0.,1.27) # Use only left side of f_2(1270)
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
