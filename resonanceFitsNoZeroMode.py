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
mPi      = 0.13957018
mK       = 0.493677

mRho     =  .77549
Grho     =  .1491

def getRhoModel(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap, phaseFile, Jxi, L):
	rhoMass  = ptc.parameter( mRho,  "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, Jxi, L, False)
	fitRho = amplitudeAnalysis(inFileName, [sector], {sector:[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(phaseFile = phaseFile)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.mode = AMPL
	return fitRho
	

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
	if not sect in ['1++1+', '2-+1+', '2++1+']:
		raise RuntimeError("Invalid sector '" + sect + "'")

	if len(sys.argv) > 3:
		study = sys.argv[2]
		studyAdder = "_"+study
	else:
		study = "std11"
		studyAdder = ""
	print "Study: "+study


#	sect = '1++1+'
#	sect = '2-+1+'
#	sect = '2++1+'

	if sect == '1++1+':
		sector = "1++1+[pi,pi]1--PiS"
		L      = 0
		Jxi    = 1
	if sect == '2-+1+':
		sector = "2-+1+[pi,pi]1--PiP"
		L      = 1
		Jxi    = 1
	if sect == '2++1+':
		sector = "2++1+[pi,pi]1--PiD"
		L      = 2
		Jxi    = 1



	inFileName    = fileNameMap[study]
#	phaseFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/ampls_4++1+rhoPiG.dat"
	phaseFileName = ""

	tBins            = [tBin]
	startBin         = 11
	stopBin          = 50

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

	model = getRhoModel(inFileName, sector, startBin, stopBin, tBins, {}, phaseFile = phaseFileName, Jxi = Jxi, L = L)

##### Writing starts here

	parameterDummy = [[]*(stopBin-startBin)]

	doRhoFits = True
	if doRhoFits:
		with open("./rhoMassesAndWidths_rhoRange0500-1100" + sect + "_"+str(tBin)+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				while True:
					try:
#					if True:
						x,err = model.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0],[i], zeroModeParameters = parameterDummy)
						break
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							raise Exception("Too many failed attempts: "+str(exceptCount))
	
				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				outFile.write('\n')			
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
