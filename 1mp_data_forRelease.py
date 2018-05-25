import matplotlib
matplotlib.use("Agg")

from analysisClass import amplitudeAnalysis
import parameterTrackingParameterizations as ptc
import parameterizationClasses as pc
from fixedparameterizationPaths import getFileNameForSector
from modes import PHASE, AMPL, SMOOTH, NONE
from utils import sumUp, weightedSum, cloneZeros, checkLaTeX
import sys, os
import pyRootPwa
import numpy as np
from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0
from rootfabi import root_open
import LaTeX_strings

import consistencyUtils as cu

import modernplotting.colors
import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from LaTeX_strings import unCorrected_string, weightedAVG_string

def doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None, avc = None):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"1-+1+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.removeZeroModeFromComa()
	if acv is not None:
		fitRho.addComaValueForZeroMode(acv)
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
	if writeResultToFile is not None:
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

	tBins            = [tBin]

#	startBin         = 11
#	stopBin          = 50

	startBin = 26
	stopBin  = 28

#	acv = 130. # artificial coma value
	acv = 1.3e-5 # artificial coma value

	methodBinRanges = {
	                   "fitF2"        : (22, 50),
	                   "fitRhoF"      : (22, 50),
	                   "fixedShapeF2" : (22, 50)}
#	methodBinRanges = {} # Override here 


	sectorRangeMap = {"1-+1+[pi,pi]1--PiP":(0.,1.12)}	
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


	scalle     = False
	Xcheck     = True
	resolvedWA = fixedShapes.getZeroModeParametersForMode()
	if scalle:
		folder = "./comparisonResultsData_1mp_scale/"
	else:
		folder = "./comparisonResultsData_1mp/"

#	folder = "./ellipseComparisons/"
	folder = "./"

	fileSectorName = "1mp1p1mmP"

	for s, sect in enumerate(fullRanges.sectors):
		fullRanges.removeZeroModeFromComa()
#		rv = fullRanges.produceResultViewer(cloneZeros(resolvedWA),s, noRun = True, plotTheory = True)
		rv = fullRanges.produceResultViewer(resolvedWA,s, noRun = True, plotTheory = True)
#		rv.writeToRootFile("exoticForRealease_"+str(tBin)+".root")

		rv.titleFontSize = 11
		rv.showColorBar  = True
		rv.XcheckArgand  = Xcheck
		rv.writeBinToPdf(startBin, stdCmd = [ fileSectorName + "_2D_t"+str(tBin)+".pdf", "", [], "", []])

		rv.labelPoints     = [0,10,15,20]
		if tBin in [0,1,2] and not scalle:
			rv.labelPoints     = [10,15,20]
		rv.makeLegend      = True
		if scalle:
			rv.scaleTo = "maxCorrData"
			rv.yAxisShift      = 300.
			rv.tStringYpos     = 0.8
			rv.topMarginIntens = 1.4
			fakkkk             = 1.
		else:
			rv.plotData        = False
			fakkkk             = .7	
			rv.tStringYpos     = 0.8
			rv.topMarginIntens = 1.4
			rv.yAxisShift      = 100.
			rv.scaleArgandZeroLine = 0.8

		rv.addiColor     = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, .5)
		rv.realLabel     = LaTeX_strings.realReleaseNote
		rv.imagLabel     = LaTeX_strings.imagReleaseNote
		rv.m2PiString    = LaTeX_strings.m2PiReleaseNote
		rv.intensLabel   = LaTeX_strings.intensReleaseNote
		rv.printLiminary = True
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

		rv.scale(1000./40.)
		for b in range(startBin, stopBin):
			if not b == 27:
				continue
#			intensNames = [name+".intens" for name in fileNames[sect,b]]
#			argandNames = [name+".argand" for name in fileNames[sect,b]]
			if not scalle:
				intensNames = []
				argandNames = []
#			intensNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/1mp1p1mm_t"+str(tBin)+".intens"]
#			argandNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/1mp1p1mm_t"+str(tBin)+".argand"]

			intensNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/filesForRelease/fullRangeFixing_"+str(tBin)+".intens"]
			argandNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/filesForRelease/fullRangeFixing_"+str(tBin)+".argand"]

			if Xcheck:
				intensNames    = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/1mp1p1mm_t"+str(tBin)+".intens"]
				argandBaseName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/dimasEllipses/1mp1p1mm_"+str(b)+"_"+str(tBin)+"_<ndx>.argand"
				argandNames    = [argandBaseName.replace("<ndx>", ndx) for ndx in ['0','1','2','3','C']]

#			print argandNames
#			return

#			intensNames = []
#			argandNames = []

			for fn in intensNames:
				if not os.path.isfile(fn):
					raise IOError
			for fn in argandNames:
				if not os.path.isfile(fn):
					raise IOError

#			rv.writeAmplFiles(b, fileName = "./filesForRelease/fullRangeFixing_"+str(tBin))

			rv.legendMethods    = "Full range"
			if Xcheck:
				rv.addiXoffset      = 0.01
				rv.addiColor        = modernplotting.colors.colorScheme.red
				rv.legendMethods    = "X-check"

#			rv.tString          = ""
#			rv.scaleFakk        = 1000.
#			rv.writeBinToPdf(b, stdCmd = ["","",[], folder+  sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
			rv.writeBinToPdf(b, stdCmd = ["",folder +  fileSectorName + "_int_m"+str(b)+"_t"+str(tBin)+".pdf", intensNames, folder+  fileSectorName + "_arg_m"+str(b)+"_t"+str(tBin)+".pdf", argandNames])
#			rv.wiriteReImToPdf(b, sect + "_data_<ri>_"+str(b)+"_"+str(tBin)+".pdf")

	with root_open("1mp1p_totals_t"+str(tBin)+".root", "RECREATE"):
		hists = fullRanges.makeTheoryTotals()
		for tb in hists:
			for h in tb:
				h.Write()
		totalHists = fixedShapes.getTotalHists(resolvedWA)
		for tb in totalHists:
			for h in tb:
				h.Write()
	return

if __name__ == "__main__":
	main()
