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

import modernplotting.colors
import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter
from globalDefinitions import mPi, mK, mRho, Grho, mRhoPrime, GrhoPrime, mF0, g1, g2, mF2, GF2, Pr0
from globalDefinitions import referenceWave as globalReferenceWave

def getRhoModel(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap,  L, referenceWave = "", writeResultToFile = None, loadIntegrals = False):
	rhoMass  = ptc.parameter( mRho-.1,  "rhoMass" )
	rhoWidth = ptc.parameter( Grho+.1 , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, L, False)
	fitRho = amplitudeAnalysis(inFileName, [sector], {sector:[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave, loadIntegrals = loadIntegrals)
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


def getF2Model(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap,  L, referenceWave = "", writeResultToFile = None, loadIntegrals = False):
	print 'sector', sector,'has L =',L
        f2Mass  = ptc.parameter( mF2,  "f2Mass" )
        f2Width = ptc.parameter( GF2 , "f2Width")
        f2 = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, L, False)
        fitF2 = amplitudeAnalysis(inFileName, [sector], {sector:[f2]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
        fitF2.loadData(referenceWave = referenceWave, loadIntegrals = loadIntegrals)
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

	referenceWave = globalReferenceWave

	tBin = int(sys.argv[1])

	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	sect = sys.argv[2]
	if not sect in ['1++1+', '2-+1+','2-+1+2++', '2-+1+2++', '2++1+','2++1+2++', '4++1+1--', '4++1+2++', '3++0+', "4-+0+", "6-+0+"]:
		raise RuntimeError("Invalid sector '" + sect + "'")

	if len(sys.argv) > 3:
		study = sys.argv[3]
		studyAdder = "_"+study
	else:
		study = "std11"
		studyAdder = ""
	print "Study: "+study

	intensNames            = []
	argandNames            = []
	Xcheck                 = False

#	sect                   = '1++1+'
#	sect                   = '2-+1+'
#	sect                   = '2++1+'
	startBin               = 11
	stopBin                = 50
	isRho                  = True
	sectorRangeMap         = {}
	isobName               = "rho"
	writBins               = None
	labelPoints            = {0:[],1:[],2:[],3:[]}
	shiftMap               = {0:{},
	                          1:{},
	                          2:{},
	                          3:{}}
	tickMap                = {}
	intensFakk             = {0:1., 1:1., 2:1., 3:1.}
	if sect == '1++1+':
		sector         = "1++1+[pi,pi]1--PiS"
		L              = 0
	if sect == '2-+1+':
		sector         = "2-+1+[pi,pi]1--PiP"
		L              = 1
	if sect == '2-+1+2++':
		sector         = "2-+1+[pi,pi]2++PiS"
		L              = 0
		isRho          = False
#		startBin       = 22
		isobName       = "f2"
	if sect == '2++1+':
		sector         = "2++1+[pi,pi]1--PiD"
		L              = 2
	if sect == "2++1+2++":
		sector         = "2++1+[pi,pi]2++PiP"
		L              = 1
		stRtBin        = 22
		isRho          = False
		isobName       = "f2"
	if sect == '4++1+1--':
		sector         = '4++1+[pi,pi]1--PiG'
		L              = 4
		writeBins      = [36]
		labelPoints    = {0:[10, 15, 20], 1:[10, 15, 21], 2:[10, 15, 20], 3:[10, 15, 20]}
		shiftMap       = {0:{10:(5*20.,0.), 15: (5*-50., 0.), 20: (0., 5*-40.)},
		                  1:{10:(0., 2.*1000.**.5), 15:(0.1*1000., 0.), 21: (0., -150.)},
		                  2:{10:(55.,0.),15:(-220.,50.),20:(0.,50.)},
		                  3:{10:(0.,50.),15:(80.,0.),20:(0.,50.)}}
		intensFakk     = {0:1.05, 1:1.1, 2:1., 3:1.1}
#		sectorRangeMap = {'4++1+[pi,pi]1--PiG':(0., 1.2)}
		if Xcheck:
			intensNames    = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/4pp1p1mm_t0.intens"]
			argandNames    = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/4pp1p1mm_t0.argand"]
	if sect == '4++1+2++':
		sector         = '4++1+[pi,pi]2++PiF'
		L              = 3
		isRho          = False	
		writeBins      = [36]
		labelPoints    = {0:[25, 30], 1:[ 25, 30], 2:[25, 30], 3:[25, 30]}
		shiftMap       = {0:{25:(50.,0.), 30:(-130.,0.)},
		                  1:{25:(0.,25.), 30:(-150., 0.)},
		                  2:{25:(0.,-60.), 30:(0., 60.)},
		                  3:{25:(0.,-80.), 30:(40., 0.)}}
		tickMap[1]     = ([-400, -200, 0, 200], [0., 200., 400., 600.])
#		startBin       = 22
		isobName       = "f2"
		if Xcheck:
			intensNames    = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/4pp1p2pp_t0.intens"]
			argandNames    = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/singlePlots/4pp1p2pp_t0.argand"]
	if sect == "3++0+":
		sector         = "3++0+[pi,pi]2++PiP"
		L              = 1
		isRho          = False
		writeBins      = [34]
		labelPoints    = {0:[25, 30], 1:[ 25, 30], 2:[25, 30], 3:[25, 30]}
		shiftMap       = {0:{25:(15.,15.), 30:( 50.,-40.)},
		                  1:{25:(10.,10.), 30:( 10., 10.)},
		                  2:{25:(15.,15.), 30:(-90., 40.)},
		                  3:{25:(10.,10.), 30:( 10., 10.)}}
#		startBin       = 22
		intensFakk     = {0:1.1, 1:1.1, 2:1.1, 3:1.2}
		isobName       = "f2"
	if sect == "4-+0+":
		sector         = "4-+0+[pi,pi]1--PiF"
		referenceWave  = "6-+0+rhoPiH"
		L              = 3
	if sect == "6-+0+":
		sector         = "6-+0+[pi,pi]1--PiH"
		L              = 5



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
		model = getRhoModel(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap, L = L, writeResultToFile = "rhoMassesAndWidths_"+sect+"_global.dat", referenceWave = referenceWave, loadIntegrals = True)
	else:
		print "Fit is f2"
		model = getF2Model(inFileName, sector, startBin, stopBin, tBins, sectorRangeMap, L = L, writeResultToFile = "f2MassesAndWidths_"+sect+"_global.dat", referenceWave = referenceWave, loadIntegrals = True)
##### Writing starts here
	print "From",startBin,'to',stopBin

	parameterDummy = [[]*(stopBin-startBin)]
	if isRho:
		outFileName = "./rhoMassesAndWidths"
	else:
		outFileName = "./f2MassesAndWidths"
	doResonanceFits = False
	binWiseFit      = True
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
					try:
#					if True:
						if binWiseFit:
							fitRange = [i]
						else:
							fitRange = range(stopBin-startBin)
						if isRho:
							startVals = [mRho+startValueOffset,Grho+startValueOffset]
						else:
							startVals = [mF2+startValueOffset,GF2+startValueOffset]
						x,err, c2, ndf = model.fitShapeParametersForBinRange(startVals, [0],fitRange, zeroModeParameters = parameterDummy)
						break
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							raise Exception("Too many failed attempts: "+str(exceptCount))
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

	makeTotals = False
	if makeTotals:
		model.calculateNonShapeParameters()
		totalHists = model.getTotalHists(model.getZeroModeParametersForMode())
	
		with root_open("./totals_2mp1p"+studyAdder+".root", "UPDATE") as out:
			for t in totalHists:
				for m in t:
					m.Write()

	folder = "./comparisonResultsData"+studyAdder+"/"
	for s, sect in enumerate(model.sectors):
		model.removeZeroModeFromComa()
		model.removeGlobalPhaseFromComa()
		zeroParDummy = [[] * (stopBin-startBin)]
		model.calculateNonShapeParameters()
		rv = model.produceResultViewer(zeroParDummy,s, noRun = True, plotTheory = True)
		rv.plotData = False
		rv.writeBinToPdf(startBin, stdCmd = [folder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
		if not writeBins:
			writeBins = range(startBin, stopBin)

		for b in writeBins:
#			intensNames         = [name+".intens" for name in fileNames[sect,b]]
#			argandNames         = [name+".argand" for name in fileNames[sect,b]]
			rv.plotData         = False
			rv.realLabel        = LaTeX_strings.realReleaseNote
			rv.imagLabel        = LaTeX_strings.imagReleaseNote
			rv.m2PiString       = LaTeX_strings.m2PiReleaseNote
			rv.tStringYpos      = 0.865
			rv.addiColor        = modernplotting.colors.colorScheme.red
			rv.legendCorrected  = r"Freed-isobar PWA"
			if Xcheck:			
				rv.legendMethods    = "X-check"
				rv.tString          = ""
			rv.labelPoints      = labelPoints[tBin]
			rv.shiftMap         = shiftMap[tBin]
			rv.makeLegend       = True
			rv.topMarginIntens *= intensFakk[tBin]
			rv.intensLabel      = LaTeX_strings.intensReleaseNote
			rv.printLiminary    = True
			rv.scale(1000./40.)
			if tBin in tickMap:
				rv.xticks = tickMap[tBin][0]
				rv.yticks = tickMap[tBin][1]

			rv.writeBinToPdf(b, stdCmd = ["", folder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  folder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])


if __name__ == "__main__":
	main()
