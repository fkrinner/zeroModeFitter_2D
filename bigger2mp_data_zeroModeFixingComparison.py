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

import consistencyUtils as cu
import LaTeX_strings

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter
mPi      = 0.13957018
mK       = 0.493677

mRho     =  .77549
Grho     =  .1491

mF2 = 1.2751
GF2 = 0.1851

def doFitRhoP(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitRho.loadData()
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitRhoF(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiF":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitRho.loadData()
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitBothRho(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rhoP = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	rhoF = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiP":[rhoP], "2-+0+[pi,pi]1--PiF" : [rhoF]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitRho.loadData()
	fitRho.finishModelSetup()
#	for mBin in fitRho.model[0]:
#		print ">>",mBin.nPar,"<<"
#	raise Exception
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitBothF2(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	f2Mass  = ptc.parameter( mF2,  "f2Mass" )
	f2Width = ptc.parameter( GF2 , "f2Width")
	f2S = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, False)
	f2D = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 2, False)
	fitF2 = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]2++PiS":[f2S], "2-+0+[pi,pi]2++PiD":[f2D]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitF2.loadData()
	fitF2.finishModelSetup()
	fitF2.fitShapeParameters()
	fitF2.calculateNonShapeParameters()
	fitF2.mode = AMPL
#	fitF2.removeGlobalPhaseFromComa()
	return fitF2



alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def doFixedShapes(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
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
	fixedShapes.loadData(loadIntegrals = True)
	fixedShapes.finishModelSetup()
	fixedShapes.fitShapeParameters()
	fixedShapes.calculateNonShapeParameters()
	fixedShapes.mode = AMPL
#	fixedShapes.removeGlobalPhaseFromComa()
	return fixedShapes

def doSmooth(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	waveModel = {}
	smooth = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	smooth.loadData()
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
	inFileName   = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_bigger2mp.root"
	zeroFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/zeroModes_bigger2mp.root"
	sectors          = ["2-+0+[pi,pi]0++PiD", "2-+0+[pi,pi]1--PiP", "2-+0+[pi,pi]1--PiF", "2-+0+[pi,pi]2++PiS","2-+0+[pi,pi]2++PiD"]
	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	tBins            = [tBin]

	startBin         = 11
	stopBin          = 50

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
	                      "smooth"            : r"smooth"}

	print "Starting with fixed shape f0"
	fixedShapeF0 = doFixedShapes(inFileName,zeroFileName, sectors[:1], startBin, stopBin, tBins)
	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"

#	print "Starting with fixed shape rhoP"
#	fixedShapeRhoP = doFixedShapes(inFileName, [sectors[1]], startBin, stopBin, tBins)
#	allMethods["fixedShapeRhoP"] = fixedShapeRhoP
#	print "Finished with fixed shape rhoP"

#	print "Starting with fixed shape rhoF"
#	fixedShapeRhoF = doFixedShapes(inFileName, [sectors[2]], startBin, stopBin, tBins)
#	allMethods["fixedShapeRhoF"] = fixedShapeRhoF
#	print "Finished with fixed shape rhoF"

#	print "Starting with fixed shape bothRho"
#	fixedShapeBothRho = doFixedShapes(inFileName, sectors[1:3], startBin, stopBin, tBins)
#	allMethods["fixedShapeBothRho"] = fixedShapeBothRho
#	print "Finished with fixed shape bothRho"

#	print "Starting with fixed shape f2"
#	fixedShapeF2 = doFixedShapes(inFileName, [sectors[3]], startBin, stopBin, tBins)
#	allMethods["fixedShapeF2"] = fixedShapeF2
#	print "Finished with fixed shape f2"

	print "Starting with fixed shapes"
	fixedShapes = doFixedShapes(inFileName,zeroFileName, sectors, startBin, stopBin, tBins)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"
	fullSig = fixedShapes.getZeroModeSignature()

#	totalHists = fixedShapes.getTotalHists(cloneZeros(fixedShapes.getZeroModeParametersForMode()))
#	with root_open("./totals_noCorr_bigger1pp.root", "UPDATE") as out:
#		for t in totalHists:
#			for m in t:
#				m.Write()
#	return

##	print "Starting with phase"
##	fitPiPiSshape = doF0phase(inFileName, sectors[:1], startBin, stopBin, tBins)
##	allMethods["pipiS"] = fitPiPiSshape
##	print "Finished with phase"

#	print "Starting with fitting rhoP"
#	fitRhoP = doFitRhoP(inFileName, sectors, startBin, stopBin, tBins)
#	allMethods["fitRhoP"] = fitRhoP
#	print "Finished with fitting rhoP"

#	print "Starting with fitting rhoF"
#	fitRhoF = doFitRhoF(inFileName, sectors, startBin, stopBin, tBins)
#	allMethods["fitRhoF"] = fitRhoF
#	print "Finished with fitting rhoF"

	print "Starting with fitting bothRho"
	fitBothRho = doFitBothRho(inFileName,zeroFileName, sectors, startBin, stopBin, tBins)
	allMethods["fitBothRho"] = fitBothRho
	print "Finished with fitting bothRho"

	print "Starting with fitting bothf2"
	fitBothF2 = doFitBothF2(inFileName,zeroFileName, sectors, startBin, stopBin, tBins)
	allMethods["fitBothF2"] = fitBothF2
	print "Finished with fitting bothf2"

	if stopBin - startBin > 1:
		print "Starting with smooth"
		smooth = doSmooth(inFileName,zeroFileName, sectors, startBin, stopBin, tBins)
		allMethods["smooth"] = smooth
		print "Finished with smooth"

	if "fixedShapeF0" in allMethods:
		allMethods["fixedShapeF0"].setZeroModeSignature(fullSig,1)

	diffsFull, resolvedWA, nonResolvedWA, comps, resDiffs, nonResDiffs, resolvedDiffsFull,noCorrDiffs = cu.doAllComparisons(allMethods, startBin, methodBinRanges)
#	print resolvedDiffsFull
	from math import isnan
	for pair in resolvedDiffsFull:
		with  modernplotting.toolkit.PdfWriter('./resolvedDiffPlots/bigger2mp_'+pair[0]+"_"+pair[1]+"_"+str(tBin)+".pdf") as pdfOutput:
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

	with  modernplotting.toolkit.PdfWriter("compositions_bigger2mp_"+str(tBin)+".pdf") as pdfOutput:
		plot = style.getPlot1D()
		for m in studyList:
			line  = [0.]*len(comps[m][0])
			xAxis = [ .5 + 0.04*(startBin + i) for i in range(len(comps[m][0]))]
			break
		count = 0
		for m in studyList:
			newLine = line[:]
			for i in range(len(comps[m][0])):
				newLine[i] += comps[m][0][i]
			plot.axes.fill_between(xAxis, line, newLine, facecolor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, 0.1*count))
			count += 1
			line = newLine
		plot.setYlim(0.,1.)
		plot.setXlim(xAxis[0], xAxis[-1])
		pdfOutput.savefigAndClose()



	hist = pyRootPwa.ROOT.TH2D("hist","hist", len(studyList)+2, 0, len(studyList)+2, len(studyList), 0, len(studyList))

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

	with modernplotting.toolkit.PdfWriter("studies_bigger2mp_"+str(tBin)+".pdf") as pdfOutput:
		plot = style.getPlot2D()
		plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(len(studyList)+2)])
		plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(len(studyList))])
		studyPlotter.makeValuePlot(plot, hist)
		
		plot.axes.set_yticklabels(axolotl)
		axolotl.append(r"$\vec 0$")
		axolotl.append(r"$\Omega$")
		plot.axes.set_xticklabels(axolotl, rotation = 90)
		plot.setZlim((0.,1.))

		pdfOutput.savefigAndClose()

	with open("studies_bigger2mp_"+str(tBin)+".txt",'w') as out:
		for axl in axolotl:
			out.write(axl + ' ')
		out.write("\n")
		for i in range(hist.GetNbinsX()):
			for j in range(hist.GetNbinsY()):
				out.write(str(hist.GetBinContent(i+1, j+1)) + ' ')
			out.write('\n')


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
				fileName = "./collectedMethods/"+stu+"_"+sect+"_bigger2mp_"+str(bin)
				if not (sect,bin) in fileNames:
					fileNames[sect,bin] = []
				fileNames[sect,bin].append(fileName)
				rv.writeAmplFiles(bin, fileName = fileName)

	totalHists = fixedShapes.getTotalHists(resolvedWA)
	with root_open("./totals_bigger2mp.root", "UPDATE") as out:
		for t in totalHists:
			for m in t:
				m.Write()
#	return

	plotFolder = "./comparisonResultsBigger2mp/"

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
		allMethods['fixedShapes'].removeZeroModeFromComa()
		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = allMethods['fixedShapes'].produceResultViewer(resolvedWA,s, noRun = True, plotTheory = True)
		rv.writeBinToPdf(startBin, stdCmd = [plotFolder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
		for b in range(startBin, stopBin):
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			rv.writeBinToPdf(b, stdCmd = ["", plotFolder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  plotFolder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
	print studyList


if __name__ == "__main__":
	main()
