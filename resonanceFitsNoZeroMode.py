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

from rootfabi      import root_open
import LaTeX_strings

import consistencyUtils as cu
from studyFileNames import fileNameMap

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter
from globalDefinitions import referenceWave, mPi, mK, mRho, Grho, mRhoPrime, GrhoPrime, mF0, g1, g2, mF2, GF2, Pr0

def getRhoModel(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap,  L, referenceWave = "", writeResultToFile = None):
	rhoMass  = ptc.parameter( mRho,  "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, L, False)
	fitRho = amplitudeAnalysis(inFileName, [sector], {sector:[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.mode = AMPL
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(rhoMass.value) + ' ' + str(rhoMass.error) + ' ' + str(rhoWidth.value) + ' ' + str(rhoWidth.error) + "\n"
			outFile.write(resultString)		
	return fitRho


def getF2Model(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap,  L, referenceWave = "", writeResultToFile = None):
	print 'sector', sector,'has L =',L
        f2Mass  = ptc.parameter( mF2,  "f2Mass" )
        f2Width = ptc.parameter( GF2 , "f2Width")
        f2 = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, L, False)
        fitF2 = amplitudeAnalysis(inFileName, [sector], {sector:[f2]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
        fitF2.loadData(referenceWave = referenceWave)
        fitF2.finishModelSetup()
        fitF2.fitShapeParameters()
        fitF2.mode = AMPL
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(f2Mass.value) + ' ' + str(f2Mass.error) + ' ' + str(f2Width.value) + ' ' + str(f2Width.error) + "\n"
			outFile.write(resultString)
        return fitF2
	

def main():
	checkLaTeX()
	style = modernplotting.mpplot.PlotterStyle()
#	style.p2dColorMap = 'ocean_r'
#	style.p2dColorMap = 'YlOrRd'
	style.p2dColorMap = 'Reds'

	tBin = int(sys.argv[1])

	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	sect = sys.argv[2]
	if not sect in ['1++1+', '2-+1+', '2++1+', '4++1+1--', '4++1+2++', '3++0+']:
		raise RuntimeError("Invalid sector '" + sect + "'")

	if len(sys.argv) > 3:
		study = sys.argv[3]
		studyAdder = "_"+study
	else:
		study = "std11"
		studyAdder = ""
	print "Study: "+study


#	sect = '1++1+'
#	sect = '2-+1+'
#	sect = '2++1+'
	startBin         = 11
	stopBin          = 50
	isRho            = True
	sectorRangeMap   = {}
	isobName         = "rho"
	if sect == '1++1+':
		sector   = "1++1+[pi,pi]1--PiS"
		L        = 0
	if sect == '2-+1+':
		sector   = "2-+1+[pi,pi]1--PiP"
		L        = 1
	if sect == '2++1+':
		sector   = "2++1+[pi,pi]1--PiD"
		L        = 2
	if sect == '4++1+1--':
		sector   = '4++1+[pi,pi]1--PiG'
		L        = 4
#		sectorRangeMap = {'4++1+[pi,pi]1--PiG':(0., 1.2)}
	if sect == '4++1+2++':
		sector   = '4++1+[pi,pi]2++PiF'
		L        = 3
		isRho    = False	
		startBin = 22
		isobName = "f2"
	if sect == "3++0+":
		sector   = "3++0+[pi,pi]2++PiP"
		L        = 1
		isRho    = False
		startBin = 22
		isobName = "f2"



	inFileName    = fileNameMap[study]
#	phaseFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/ampls_4++1+rhoPiG.dat"
	phaseFileName = ""

	tBins            = [tBin]

	allMethods       = {}
	methodStrings    = {}
	shortlabels      = {  "fixedShapeF0"    : r"$\text{fix}_{f_0}^{~}$",
	                      "fixedShapeRho"   : r"$\text{fix}_\rho^{~}$",
	                      "fixedShapeRho1G" : r"$\text{fix}_\rho^{1\Gamma}$",
	                      "fixedShapeRho2G" : r"$\text{fix}_\rho^{2\Gamma}$",
	                      "fixedShapes"     : r"$\text{fix}_\text{all}^{~}$",
	                      "pipiS"           : r"$\phi_{[\pi\pi]_S}^{~}$",
	                      "fitRho"          : r"$\text{fit}_\rho^{~}$",
	                      "fitRho1G"        : r"$\text{fit}_\rho^{1\Gamma}$",
	                      "fitRho2G"        : r"$\text{fit}_\rho^{2\Gamma}$",
	                      "smooth"          : r"smooth"}

	if isRho:
		print "Fit is rho"
		model = getRhoModel(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap, L = L, writeResultToFile = "rhoMassesAndWidths_"+sect+"_global.dat", referenceWave = referenceWave)
	else:
		print "Fit is f2"
		model = getF2Model(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap, L = L, writeResultToFile = "f2MassesAndWidths_"+sect+"_global.dat", referenceWave = referenceWave)
##### Writing starts here
	print "From",startBin,'to',stopBin

	parameterDummy = [[]*(stopBin-startBin)]
	if isRho:
		outFileName = "./rhoMassesAndWidths"
	else:
		outFileName = "./f2MassesAndWidths"
	doResonanceFits = True
	binWiseFit      = False
	if doResonanceFits:
		with open(outFileName + '_' + sect + "_"+str(tBin)+".dat",'w') as outFile:
			if binWiseFit:
				bins = range(stopBin-startBin)
			else:
				bins = [1]
			for i in bins:
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				while True:
#					try:
					if True:
						if binWiseFit:
							fitRange = [i]
						else:
							fitRange = range(stopBin-startBin)
						x,err, c2, ndf = model.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0],fitRange, zeroModeParameters = parameterDummy)
						break
#					except:
#						print "Fitter exception encountered"
#						startValueOffset += 0.001
#						exceptCount      += 1	
#						if exceptCount > 3:
#							raise Exception("Too many failed attempts: "+str(exceptCount))
				if sect == "3++0+":
					with open("3pp_f2_cpls_"+str(tBin)+".dat",'w') as outFileCpl:
						model.calculateNonShapeParameters()
						cpl, hess = model.getCouplingParameters()
						outFileCpl.write(str(cpl))
						outFileCpl.write("\n")
						outFileCpl.write(str(hess))
				nBin = startBin + i
				if nBin in []:
					model.calculateNonShapeParameters()
					folder = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/zeroModeFitter_2D/rhoShapes4pp"
					rv = model.produceResultViewer(parameterDummy,sector, noRun = True, plotTheory = True)
					rv.writeBinToPdf(nBin, stdCmd = ["",folder +  "/"+isobName+"Shape_"+sector+"_intens_"+str(nBin)+"_"+str(tBin)+".pdf", [], folder +  "/"+isobName+"Shape_"+sector+"_argand_"+str(nBin)+"_"+str(tBin)+".pdf", []])
				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				outFile.write(' ' + str(c2/ndf) + '\n')			
		return

	fileNames = {}

	for stu in allMethods:
		print "Writing for '" + stu + "'"
		for s, sect in enumerate(allMethods[stu].sectors):
			if stu == "pipiS":
				rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, plotTheory = True)
				rv.run()
			rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, noRun = True)
			for bin in range(startBin, stopBin):
				fileName = "./collectedMethods/"+stu+"_"+sect+"_0mpData_"+str(bin)+studyAdder
				if not (sect,bin) in fileNames:
					fileNames[sect,bin] = []
				fileNames[sect,bin].append(fileName)
				rv.writeAmplFiles(bin, fileName = fileName)

	totalHists = fixedShapes.getTotalHists(resolvedWeightedSum)
	
	with root_open("./totals_0mp"+studyAdder+".root", "UPDATE") as out:
		for t in totalHists:
			for m in t:
				m.Write()
#	return
	folder = "./comparisonResultsData"+studyAdder+"/"
	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
		allMethods['fixedShapes'].removeZeroModeFromComa()
		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = allMethods['fixedShapes'].produceResultViewer(resolvedWeightedSum,s, noRun = True, plotTheory = True)
		rv.writeBinToPdf(startBin, stdCmd = [folder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
		for b in range(startBin, stopBin):
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			rv.writeBinToPdf(b, stdCmd = ["", folder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  folder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
	print studyList
	print cumulWeights


if __name__ == "__main__":
	main()
