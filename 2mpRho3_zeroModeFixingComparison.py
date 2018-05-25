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

from rootfabi import root_open
from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0,mRho3,Grho3
import consistencyUtils as cu
import LaTeX_strings

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from LaTeX_strings import unCorrected_string, weightedAVG_string

def doFitRhoP(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitRhoF(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiF":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitBothRho(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rhoP = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	rhoF = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiP":[rhoP], "2-+0+[pi,pi]1--PiF" : [rhoF]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
#	for mBin in fitRho.model[0]:
#		print ">>",mBin.nPar,"<<"
#	raise Exception
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitBothF2(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	f2Mass  = ptc.parameter( mF2,  "f2Mass" )
	f2Width = ptc.parameter( GF2 , "f2Width")
	f2S = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, False)
	f2D = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 2, False)
	fitF2 = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]2++PiS":[f2S], "2-+0+[pi,pi]2++PiD":[f2D]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitF2.loadData(referenceWave = referenceWave)
	fitF2.finishModelSetup()
	fitF2.fitShapeParameters()
	fitF2.calculateNonShapeParameters()
	fitF2.mode = AMPL
#	fitF2.removeGlobalPhaseFromComa()
	return fitF2

def doFitF2S(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None):
	f2Mass  = ptc.parameter( mF2, "f2Mass" )
	f2Width = ptc.parameter( GF2 , "f2Width")
	f2 = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]2++PiS":[f2]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
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
			resultString = str(tBins[0])+ " 666. " + str(f2Mass.value) + ' ' + str(f2Mass.error) + ' ' + str(f2Width.value) + ' ' + str(f2Width.error) + "\n"
			outFile.write(resultString)
	return fitRho

def doFitF2D(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None):
	f2Mass  = ptc.parameter( mF2, "f2Mass" )
	f2Width = ptc.parameter( GF2 , "f2Width")
	f2 = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 2, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]2++PiD":[f2]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
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
			resultString = str(tBins[0])+ " 666. " + str(f2Mass.value) + ' ' + str(f2Mass.error) + ' ' + str(f2Width.value) + ' ' + str(f2Width.error) + "\n"
			outFile.write(resultString)
	return fitRho

def doFitRho3(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None):
	rho3mass  = ptc.parameter( mRho3,  "rho3mass" )
	rho3width = ptc.parameter( Grho3 , "rho3width")
	rho3 = ptc.relativisticBreitWigner([rho3mass,rho3width], mPi, mPi, mPi, 3, 1, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]3--PiP":[rho3]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
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
			resultString = str(tBins[0])+ " 666. " + str(rho3mass.value) + ' ' + str(rho3mass.error) + ' ' + str(rho3width.value) + ' ' + str(rho3width.error) + "\n"
			outFile.write(resultString)
	return fitRho

alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def doFixedShapes(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	waveModel = {}
	for n,sector in enumerate(sectors):
		model = []
		fileNames = getFileNameForSector(sector, False, False)
		print fileNames
		for fn in fileNames:
			param = pc.fixedParameterization(fn, polynomialDegree  = 0, complexPolynomial = False)
			model.append(param)
		waveModel[sector] = model
	fixedShapes = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fixedShapes.loadData(loadIntegrals = True, referenceWave = referenceWave)
	fixedShapes.finishModelSetup()
	fixedShapes.fitShapeParameters()
	fixedShapes.calculateNonShapeParameters()
	fixedShapes.mode = AMPL
#	fixedShapes.removeGlobalPhaseFromComa()
	return fixedShapes

def doSmooth(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	waveModel = {}
	smooth = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
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

#	inFileName   = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_MC.root"
	inFileName   = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_rho3_2mp3pp.root"
	zeroFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_rho3_2mp3pp.root"
	sectors          = ["2-+0+[pi,pi]0++PiD", "2-+0+[pi,pi]1--PiP", "2-+0+[pi,pi]1--PiF", "2-+0+[pi,pi]2++PiS","2-+0+[pi,pi]2++PiD","2-+0+[pi,pi]3--PiP"]

	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	tBins            = [tBin]

	startBin         = 36
	stopBin          = 50

#	startBin         = 20
#	stopBin          = 25

	methodBinRanges = {
	                   "fitF2"        : (22, 50),
	                   "fitBothF2"    : (22, 50),
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
	                      "fitBothF2"         : r"$\text{fit}_{f_2}^{2}$",
	                      "fitF2"             : r"$\text{fit}_{f_2}$",
	                      "fitRho3"           : r"$\text{fit}_{\rho_3}$",
	                      "smooth"            : r"smooth"}


	print "Starting with fitting rho3"
	fitRho3 = doFitRho3(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "rho3MassesAndWidths_2mp_rho3_global.dat")
	allMethods["fitRho3"] = fitRho3
	print "Finished with fitting rho3"

	print "Starting with fixed shapes"
	fixedShapes = doFixedShapes(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"
#	fullSig = fixedShapes.getZeroModeSignature()

###	print "Starting with fitting rhoP"
###	fitRhoP = doFitRhoP(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
###	allMethods["fitRhoP"] = fitRhoP
###	print "Finished with fitting rhoP"

###	print "Starting with fitting rhoF"
###	fitRhoF = doFitRhoF(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
###	allMethods["fitRhoF"] = fitRhoF
###	print "Finished with fitting rhoF"

	addAll = True
	if addAll:
		print "Starting with fitting bothRho"
		fitBothRho = doFitBothRho(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
		allMethods["fitBothRho"] = fitBothRho
		print "Finished with fitting bothRho"

		print "Starting with fitting bothf2"
		fitBothF2 = doFitBothF2(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
		allMethods["fitBothF2"] = fitBothF2
		print "Finished with fitting bothf2"

		if stopBin - startBin > 1:
			print "Starting with smooth"
			smooth = doSmooth(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
			allMethods["smooth"] = smooth
			print "Finished with smooth"

	diffsFull, resolvedWA, nonResolvedWA, comps, resDiffs, nonResDiffs, resolvedDiffsFull, noCorrDiffs = cu.doAllComparisons(allMethods, startBin, methodBinRanges)

	makePlots = False
	if makePlots:

		plotFolder = "./comparisonResults2mp_rho3/"
		for s, sect in enumerate(sectors):
			rv = allMethods['fixedShapes'].produceResultViewer(resolvedWA,s, noRun = True, plotTheory = True, removeZM = False)
			rv.writeBinToPdf(startBin, stdCmd = [plotFolder+sect+"_data_2D_"+str(tBin)+".pdf","", [],"", []])
			rv.scaleTo = "maxCorrData"
	
			for b in range(startBin, stopBin):
				intensNames = []
				argandNames = []
				rv.writeBinToPdf(b, stdCmd = ["", plotFolder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  plotFolder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
		return

	doRho3Fits = True
	if doRho3Fits:
		with open("rho3massesAndWidths_2mpRho3_"+str(tBin)+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.01
				exceptCount      = 0
				while True:
#					try:
					if True:
						x,err,c2,ndf = fitRho3.fitShapeParametersForBinRange([mRho3+startValueOffset,Grho3+startValueOffset], [0],[i], zeroModeParameters = resolvedWA)
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
				outFile.write(' '+str(c2/ndf)+'\n')			
#				fitRho3.calculateNonShapeParametersForZeroModeParameters(resolvedWA)
#				rho3RV = fitRho3.produceResultViewer(resolvedWA,"2-+0+[pi,pi]3--PiP", noRun = True, plotTheory = True)
#				rho3RV.plotData = True
#				plotNameBase = "./f2FitPlots/2mp0p3mmPiP_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
#				rho3RV.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])

if __name__ == "__main__":
	main()
