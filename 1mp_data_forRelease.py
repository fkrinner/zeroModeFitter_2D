import matplotlib
matplotlib.use("Agg")

from analysisClass import amplitudeAnalysis
import parameterTrackingParameterizations as ptc
import parameterizationClasses as pc
from fixedparameterizationPaths import getFileNameForSector
from modes import PHASE, AMPL, SMOOTH, NONE
from utils import sumUp, weightedSum, cloneZeros, checkLaTeX
import sys, os
#import pyRootPwa
import numpy as np
from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0
from rootfabi import root_open
import LaTeX_strings

import numpy.linalg as la

import consistencyUtils as cu

import modernplotting.colors
import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from LaTeX_strings import unCorrected_string, weightedAVG_string

def parseMapFile(inFileName):
	retVal = {}
	with open(inFileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			if len(chunks) == 0:
				continue
			retVal[(int(chunks[0]), int(chunks[1]))] = [float(chunks[i]) for i in xrange(2,len(chunks))]
	return retVal

def getTotalAndNevents():
	totalFN = "./exotic_totals.txt"
	NevtsFN = "./nEvents_and_bestLL.txt"
	nEvts   = parseMapFile(NevtsFN)
	totals  = parseMapFile(totalFN)
	retVal  = {}
	for tb in totals:
		if totals[tb][0] == 0.:
			continue
		retVal[tb] = (totals[tb][0], nEvts[tb][0])
	return retVal

def loadLimits(fileName, tBin):
	retVal = {}
	with open(fileName, 'r') as inFile:
		for line in inFile.readlines():
			chunks = line.split()
			tBinIn = int(chunks[0])
			if not tBinIn == tBin:
				continue
			mBin = int(chunks[1])
			vals = [float(chunks[i]) for i in range(2, len(chunks))]
			retVal[mBin] = vals
	return retVal

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

	inFileName    = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_exotic.root"
	sectors          = ["1-+1+[pi,pi]1--PiP"]
	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))
	writeCpls     = False
	twoDmarkPlots = False
	scalle        = False
	Xcheck        = False
	globalScale   = False
	setRanges     = False
	writeROOT     = False
	rangeMin      = 0.
	rangeMax      = 1.12
	startBin      = 11
	stopBin       = 50
	acv           = 1.3e-2
	for i in range(2,len(sys.argv)):
		a = sys.argv[i]
		if a == "-scale":
			scalle = True
		elif a == "-2Dmark":
			twoDmarkPlots = True
		elif a == "-Xcheck":
			Xcheck = True
		elif a == "-cpls":
			writeCpls = True
		elif a == "-global":
			globalScale = True
		elif a.startswith("-rangeMin"):
			setRanges = True
			rangeMin  = float(a[9:])
		elif a.startswith("-rangeMax"):
			setRanges = True
			rangeMax  = float(a[9:])
		elif a.startswith("-start"):
			startBin = int(a[6:])
		elif a.startswith("-stop"):
			stopBin = int(a[5:])
		elif a == "-ROOT":
			writeROOT = True
		elif a.startswith("-acv"):
			acv = float(a[4:])
		else:
			raise RuntimeError("Could not parse option '" + a + "'")

	print "-------------------------------"
	print "tbin",tBin
	print"(startBin,stopBin) (",startBin,",",stopBin,")"
	print "writeCpls",writeCpls
	print "scale", scalle
	print "2Dmarks",twoDmarkPlots
	print "Xcheck",Xcheck
	print "globalScale",globalScale
	print "fitRange ("+str(rangeMin)+","+str(rangeMax)+")"
	print "writeROOT",writeROOT
	print "acv",acv
	print "-------------------------------"	

	tBins            = [tBin]

#	acv = None
#	acv = 130. # artificial coma value
#	acv = 1.3e-3 # artificial coma value

	methodBinRanges = {
	                   "fitF2"        : (22, 50),
	                   "fitRhoF"      : (22, 50),
	                   "fixedShapeF2" : (22, 50)}
#	methodBinRanges = {} # Override here 


	sectorRangeMap = {"1-+1+[pi,pi]1--PiP":(rangeMin, rangeMax)}	
#	sectorRangeMap = {}


	allMethods       = {}
	methodStrings    = {}
	shortlabels      = {  "fixedShapeF0"      : r"$\text{fix}_{f_0}^{~}$",
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
	print "Starting with fixed shapes"
	fixedShapes = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, sectorRangeMap = sectorRangeMap, acv = acv)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"
	fullSig = fixedShapes.getZeroModeSignature()

	fullRanges = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, sectorRangeMap = {}, acv = acv)
	zeroModeParameters = fixedShapes.getZeroModeParametersForMode()

	doRhoFits = False

	if setRanges:
		rangeString = "_rangeMin"+str(rangeMin)+"_rangeMax"+str(rangeMax)

	if writeCpls:
		print fixedShapes.nonShapeParameters
		print "Writing couplings..."
		fixedShapes.calculateNonShapeParametersForZeroModeParameters(zeroModeParameters)
		outFileCpl = open("./cplFiles/1mp_rho_cpls"+rangeString+"_"+str(tBin)+".dat",'w') 
		cpl, hess = fixedShapes.getCouplingParameters()
		if not len(cpl[0][0]) == 2:
			raise IndexError("Parameter count not 2, change implementation")
		for i in range(stopBin-startBin):
			binIndex = i+startBin
			hessInv = la.inv(hess[0][i])
			integral = fixedShapes.model[0][i].getIntegralForFunction(0, fixedShapes.model[0][i].funcs[0][0])
			outFileCpl.write(str(0.52 + binIndex*.04) + ' ' + str(cpl[0][i][0]**2 + cpl[0][i][0]**2) + ' ' + str(integral) + ' ')
			outFileCpl.write(str(cpl[0][i][0])     + ' ' + str(cpl[0][i][1])     + ' ')
			outFileCpl.write(str(hessInv[0][0]/2) + ' ' + str(hessInv[1][1]/2) + ' ' + str(hessInv[0][1]/2))
			outFileCpl.write("\n")
		outFileCpl.close()
		return

	doActuallyNotFit = False
	if doRhoFits:
		fitRho = doFitRho(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "rhoMassesAndWidths_1-+1+1--_global.dat", sectorRangeMap = sectorRangeMap, acv = acv)
		fitRho.initMinuitFunction(fitRho.getParameters())
		with open("rhoMassesAndWidths_exoticForPaper_"+str(tBin)+".dat",'w') as out:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				out.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				if doActuallyNotFit:
					print "The fit has been turned off, due to a workaround... :("
				else:
					exceptCount      = 0
#					try:
					if True:
						x,err,c2,ndf = fitRho.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0],[i], zeroModeParameters = zeroModeParameters)
#					except:
#						print "Fitter exception encountered"
#						startValueOffset += 0.001
#						exceptCount      += 1	
#						if exceptCount > 3:
#							raise Exception("Too many failed attempts: "+str(exceptCount))
	
					out.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
					out.write(' '+str(c2/ndf)+'\n')

				if writeCpls:
					fitRho.calculateNonShapeParametersForZeroModeParameters(zeroModeParameters)
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
					fitRho.calculateNonShapeParametersForZeroModeParameters(zeroModeParameters)
					rhoRV = fitRho.produceResultViewer(zeroModeParameters,"1-+1+[pi,pi]1--PiP", noRun = True, plotTheory = True)
					rhoRV.plotData = True
					plotNameBase = "./rhoFitPlots/1mp1p1mmPiP_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
					rhoRV.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
		if writeCpls:
			outFileCpl.close()

		return

	if (twoDmarkPlots and scalle) or (twoDmarkPlots and Xcheck):
		raise Exception("'twoDmarkPlots' only allowed without 'scalle' and 'Xcheck'")

	folderAddString = ""
	if scalle:
		folderAddString += "_scale"
	if Xcheck:
		folderAddString += "_Xcheck"
	if globalScale:
		folderAddString += "_global"



#	folder = "./ellipseComparisons/"
	folder = "./forPaper"+folderAddString+"/"

	if setRanges:
		folder = "./exoticFitPlots_ranges/"


	intensYlimFileName   = "./intensYlims"
	arangdLimitsFileName = "./argandLims"
	if scalle:
		intensYlimFileName   += "_scale"
		arangdLimitsFileName += "_scale"

	intensScales = loadLimits(intensYlimFileName, tBin)
	argandScales = loadLimits(arangdLimitsFileName, tBin)

	if globalScale:
		globalScales = {
					0: ( [9808549.03133*.9] ,[ -1176.89018363 , 2583.08091376 , -963.530201342 , 2764.05317273 ]), 
					1: ( [6250000.00   ] ,[ -1000. , 2311.5523738 , -950. , 2400. ]),
#					1: ( [15692071.5541] ,[ -2953.09588149 , 2311.5523738 , -1069.36767067 , 3461.97833088 ]),
					2: ( [5984192.64568*13.5/15.] ,[ -1000.84212555 , 1893.32841381 , -652.0 , 2498.12259086 ]),
					3: ( [3693705.14722*16.5/18.5] ,[ -1259.93088183 , 1277.85363449 , -404.005733504 , 1956.63757019 ])}
		for i in range(startBin, stopBin):
			intensScales[i] = globalScales[tBin][0]
			argandScales[i] = globalScales[tBin][1]

	fileSectorName = "1mp1p1mmP"
	if setRanges:
		fileSectorName += rangeString

	fullRanges.nonShapeParameters = fixedShapes.nonShapeParameters

	totalsAndNevents = getTotalAndNevents()

	for s, sect in enumerate(fullRanges.sectors):
		fullRanges.removeZeroModeFromComa()
#		rv = fullRanges.produceResultViewer(cloneZeros(zeroModeParameters),s, noRun = True, plotTheory = True)
		rv = fullRanges.produceResultViewer(zeroModeParameters,s, noRun = True, plotTheory = True)
		if writeROOT:
			rv.writeToRootFile(fileSectorName+".root")
			return 

#		rv.writeToRootFile("exoticForRealease_"+str(tBin)+".root")

		rv.checkArgandRanges = False

		rv.intYaxisXpos = -0.15
		rv.argYaxisXpos = -0.145

		rv.intensLims = intensScales
		rv.argandLims = argandScales

		rv.arrowsToTheo = [14]

		rv.titleFontSize = 11
		rv.showColorBar  = True
		rv.XcheckArgand  = Xcheck
		rv.writeBinToPdf(startBin, stdCmd = [ fileSectorName + "_2D_t"+str(tBin)+".pdf", "", [], "", []])

#		continue

		if twoDmarkPlots:
			rv.mark2DplotColor = modernplotting.colors.colorScheme.red

		rv.labelPoints     = [0,10,15,20]
		if tBin in [0,1,2] and not scalle:
			rv.labelPoints     = [10,15,20]
		rv.labelPoints     = []
		rv.makeLegend      = True
		if scalle:
			rv.scaleTo = "maxCorrData"
			rv.yAxisShift      = 300.
			rv.tStringYpos     = 0.8 # 0.8, if third legend item 
			rv.topMarginIntens = 1.4
			fakkkk             = 1.
		else:
			rv.plotData        = False
			fakkkk             = .7	
			rv.tStringYpos     = 0.87 # 0.8, if third legend item 
			rv.topMarginIntens = 1.4
			rv.yAxisShift      = 100.
#			rv.scaleArgandZeroLine = 0.8

		if twoDmarkPlots:
			rv.tStringYpos = 0.94

		if Xcheck:
			rv.tStringYpos     = 0.8 # 0.8, if third legend item 

		rv.addiColor     = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, .5)
		rv.realLabel     = LaTeX_strings.realReleaseNote
		rv.imagLabel     = LaTeX_strings.imagReleaseNote
		rv.m2PiString    = LaTeX_strings.m2PiReleaseNote
		rv.intensLabel   = LaTeX_strings.intensReleaseNote
		if not Xcheck:
			rv.printLiminary = True
		else:
			rv.connectAddGraphs = [(0,1),(2,3)]
		rv.xLims = (.2,2.2)
		if tBin == 3:
			if scalle:
				rv.shiftMap      = {0:(fakkkk*80.,fakkkk*-280.),10:(fakkkk*-360.,fakkkk*-50.), 15:(fakkkk*-390.,fakkkk*-30.), 20:(fakkkk*-500.,fakkkk*-10.)}
			else:
				rv.shiftMap      = {0:(fakkkk*100.,fakkkk*-280.),10:(fakkkk*-370.,fakkkk*-50.), 15:(fakkkk*-370.,fakkkk*-30.), 20:(fakkkk*-50.,fakkkk*90.)}
		if tBin == 2:
			if scalle:
				rv.shiftMap      = {0:(fakkkk*50.,fakkkk*-300.),10:(fakkkk*-480.,fakkkk*-50.), 15:(fakkkk*-500.,fakkkk*-30.), 20:(fakkkk*0.,fakkkk*70.)}
			else:
				rv.shiftMap      = {0:(fakkkk*50.,fakkkk*-400.),10:(fakkkk*-420.,fakkkk*-50.), 15:(fakkkk*-450.,fakkkk*-30.), 20:(fakkkk*0.,fakkkk*70.)}
		if tBin == 1:
			if scalle:
				rv.shiftMap      = {0:(fakkkk*50.,fakkkk*-400.),10:(fakkkk*-300.,fakkkk*-50.), 15:(fakkkk*-400.,fakkkk*-30.), 20:(fakkkk*60.,fakkkk*30.)}
			else:
				rv.shiftMap      = {0:(fakkkk*50.,fakkkk*-400.),10:(fakkkk*-350.,fakkkk*-50.), 15:(fakkkk*-400.,fakkkk*-30.), 20:(fakkkk*80.,fakkkk*30.)}
		if tBin == 0:
			if scalle:
				rv.shiftMap      = {0:(fakkkk*50.,fakkkk*-400.),10:(fakkkk*-450.,fakkkk*-50.), 15:(fakkkk*-450.,fakkkk*-30.), 20:(fakkkk*50.,fakkkk*30.)}
			else:
				rv.shiftMap      = {0:(fakkkk*50.,fakkkk*-400.),10:(fakkkk*-410.,fakkkk*-50.), 15:(fakkkk*-420.,fakkkk*-30.), 20:(fakkkk*70.,fakkkk*30.)}

		if not twoDmarkPlots:
			rv.scale(1000./40.)

		rv.writeBinToPdf(startBin, stdCmd = [folder + fileSectorName + "_2D_t"+str(tBin)+".pdf","", [], "",[]])
		for b in range(startBin, stopBin):
#			intensNames = [name+".intens" for name in fileNames[sect,b]]
#			argandNames = [name+".argand" for name in fileNames[sect,b]]

			print "Writing bin",b

			rv.rightString  = "{:.1f}\%".format(100*totalsAndNevents[(tBin, b)][0]/totalsAndNevents[(tBin, b)][1])
			if not scalle:
				intensNames = []
				argandNames = []
#			intensNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/1mp1p1mm_t"+str(tBin)+".intens"]
#			argandNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/1mp1p1mm_t"+str(tBin)+".argand"]

			intensNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/filesForRelease/fullRangeFixing_"+str(tBin)+".intens"]
			argandNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/filesForRelease/fullRangeFixing_"+str(tBin)+".argand"]

			intensNames = []
			argandNames = []

			if Xcheck:
				intensNames    = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/intensPlots/1mp1p1mm_m"+str(b)+"_t"+str(tBin)+".intens"]
				argandBaseName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/dimasEllipses/1mp1p1mm_"+str(b)+"_"+str(tBin)+"_<ndx>.argand"
				argandNames    = [argandBaseName.replace("<ndx>", ndx) for ndx in ['0','1','2','3','C']]

#			print argandNames
#			return

			for fn in intensNames:
				if not os.path.isfile(fn):
					raise IOError("'"+fn+"' does not exist")
			for fn in argandNames:
				if not os.path.isfile(fn):
					raise IOError("'"+fn+"' does not exist")

#			rv.writeAmplFiles(b, fileName = "./filesForRelease/fullRangeFixing_"+str(tBin))

			rv.legendMethods    = "Full range"
			rv.connectCorrPoints = range(70)
			rv.markRange = (10,15)
			if Xcheck:
				rv.addiXoffset      = 0.01
				rv.addiColor        = modernplotting.colors.colorScheme.red
				rv.legendMethods    = "X-check"

			if twoDmarkPlots:
				twoDplotName = folder +  fileSectorName + "_2D_m"+str(b)+"_t"+str(tBin)+".pdf"
				intPlotName  = ""
				argPlotName  = ""
			else:
				twoDplotName = ""
				intPlotName  = folder +  fileSectorName + "_int_m"+str(b)+"_t"+str(tBin)+".pdf"
				argPlotName  = folder+  fileSectorName + "_arg_m"+str(b)+"_t"+str(tBin)+".pdf"

				

#			rv.tString          = ""
#			rv.scaleFakk        = 1000.
#			rv.writeBinToPdf(b, stdCmd = ["","",[], folder+  sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
			rv.writeBinToPdf(b, stdCmd = [twoDplotName, intPlotName, intensNames, argPlotName, argandNames])
#			rv.wiriteReImToPdf(b, sect + "_data_<ri>_"+str(b)+"_"+str(tBin)+".pdf")
	print tBin,':','([',str(rv.InMax),'],[',str(rv.reMin),',',str(rv.reMax),',',str(rv.imMin),',',str(rv.imMax),'])'

	return

	with root_open("1mp1p_totals_t"+str(tBin)+".root", "RECREATE"):
		hists = fullRanges.makeTheoryTotals()
		for tb in hists:
			for h in tb:
				h.Write()
		totalHists = fixedShapes.getTotalHists(zeroModeParameters)
		for tb in totalHists:
			for h in tb:
				h.Write()
	return

if __name__ == "__main__":
	main()
