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
from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0
import consistencyUtils as cu

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from LaTeX_strings import unCorrected_string, weightedAVG_string

intensWeight = False
reweightInverseBW = False 

def doFitRhoP(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
#	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	binning  = [0.278, 0.32, 0.36, 0.4 , 0.44, 0.48, 0.52, 0.56, 0.6 , 0.64, 0.68, 0.7 , 0.72,
	            0.74 , 0.76, 0.78, 0.8 , 0.82, 0.84, 0.86, 0.88, 0.9 , 0.92, 0.96, 1.0 , 1.04, 
	            1.08 , 1.12, 1.16, 1.2 , 1.24, 1.28, 1.32, 1.36, 1.4 , 1.44, 1.48, 1.52, 1.56, 
	            1.6  , 1.64, 1.68, 1.72, 1.76, 1.8 , 1.84, 1.88, 1.92, 1.96, 2.0 , 2.04, 2.08, 
	            2.12 , 2.16, 2.2 , 2.24, 2.28]

	rho = ptc.integratedRelativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, binning, intensWeight = intensWeight, fitPr = False, reweightInverseBW = reweightInverseBW)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitRhoF(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
#	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, False)
	binning  = [0.278, 0.32, 0.36, 0.4 , 0.44, 0.48, 0.52, 0.56, 0.6 , 0.64, 0.68, 0.7 , 0.72,
	            0.74 , 0.76, 0.78, 0.8 , 0.82, 0.84, 0.86, 0.88, 0.9 , 0.92, 0.96, 1.0 , 1.04, 
	            1.08 , 1.12, 1.16, 1.2 , 1.24, 1.28, 1.32, 1.36, 1.4 , 1.44, 1.48, 1.52, 1.56, 
	            1.6  , 1.64, 1.68, 1.72, 1.76, 1.8 , 1.84, 1.88, 1.92, 1.96, 2.0 , 2.04, 2.08, 
	            2.12 , 2.16, 2.2 , 2.24, 2.28]

	rho = ptc.integratedRelativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, binning, intensWeight = intensWeight, fitPr = False, reweightInverseBW = reweightInverseBW)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiF":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitBothRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	#rhoP = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	#hoF = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, False)
	binning  = [0.278, 0.32, 0.36, 0.4 , 0.44, 0.48, 0.52, 0.56, 0.6 , 0.64, 0.68, 0.7 , 0.72,
	            0.74 , 0.76, 0.78, 0.8 , 0.82, 0.84, 0.86, 0.88, 0.9 , 0.92, 0.96, 1.0 , 1.04, 
	            1.08 , 1.12, 1.16, 1.2 , 1.24, 1.28, 1.32, 1.36, 1.4 , 1.44, 1.48, 1.52, 1.56, 
	            1.6  , 1.64, 1.68, 1.72, 1.76, 1.8 , 1.84, 1.88, 1.92, 1.96, 2.0 , 2.04, 2.08, 
	            2.12 , 2.16, 2.2 , 2.24, 2.28]

	rhoP = ptc.integratedRelativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, binning, intensWeight = intensWeight, fitPr = False, reweightInverseBW = reweightInverseBW)
	rhoF = ptc.integratedRelativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 3, binning, intensWeight = intensWeight, fitPr = False, reweightInverseBW = reweightInverseBW)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]1--PiP":[rhoP], "2-+0+[pi,pi]1--PiF" : [rhoF]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doFitF2(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = ""):
	f2Mass  = ptc.parameter( mF2,  "f2Mass" )
	f2Width = ptc.parameter( GF2 , "f2Width")
	binning = [ 0.278,0.32, 0.36, 0.4,  0.44, 0.48, 0.52, 0.56, 0.6,  0.64, 0.68, 
                             0.72, 0.76, 0.8,  0.84, 0.88, 0.92, 0.96, 1.0,  1.04, 1.08, 1.12, 
                             1.16, 1.18, 1.2,  1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 
                             1.38, 1.4,  1.44, 1.48, 1.52, 1.56, 1.6,  1.64, 1.68, 1.72, 1.76, 
                             1.8,  1.84, 1.88, 1.92, 1.96, 2.0,  2.04, 2.08, 2.12, 2.16, 2.2,
                             2.24, 2.28]
#	f2 = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, False)
	f2 = ptc.integratedRelativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 2, 0, binning, intensWeight = intensWeight, fitPr = False, reweightInverseBW = reweightInverseBW)
	fitF2 = amplitudeAnalysis(inFileName, sectors, {"2-+0+[pi,pi]2++PiS":[f2]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitF2.loadData(referenceWave = referenceWave)
	fitF2.finishModelSetup()
	fitF2.fitShapeParameters()
	fitF2.calculateNonShapeParameters()
	fitF2.mode = AMPL
#	fitF2.removeGlobalPhaseFromComa()
	return fitF2


def doF0phase(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"2-+0+[pi,pi]0++PiD" : (0.34, 2*mK - .04)}, referenceWave = ""):
	pipiSfileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_2mp0pSigmaPiD"
	pipiSw = pc.fixedParameterization(pipiSfileName, polynomialDegree  = 0, complexPolynomial = False)
	waveModel = {"2-+0+[pi,pi]0++PiD": [pipiSw]}
	fitPiPiSshape = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitPiPiSshape.loadData(referenceWave = referenceWave)
	fitPiPiSshape.finishModelSetup()
	fitPiPiSshape.phaseFit()
	fitPiPiSshape.mode = PHASE
#	fitPiPiSshape.removeGlobalPhaseFromComa()
#	fitPiPiSshape.unifyComa()
	return fitPiPiSshape

alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

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
	fixedShapes.loadData(referenceWave = referenceWave)
	fixedShapes.finishModelSetup()
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

	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_MC_corrected.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_std11.root"
	sectors          = ["2-+0+[pi,pi]0++PiD", "2-+0+[pi,pi]1--PiP", "2-+0+[pi,pi]1--PiF", "2-+0+[pi,pi]2++PiS"]
	tBins            = [0]
	startBin         = 11
	stopBin          = 50

	methodBinRanges = {
	                   "fitF2"        : (22, 50),
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
	                      "fitF2"             : r"$\text{fit}_{f_2}$",
	                      "smooth"            : r"smooth"}

	print "Starting with fixed shape f0"
	fixedShapeF0 = doFixedShapes(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"

#	print "Starting with fixed shape rhoP"
#	fixedShapeRhoP = doFixedShapes(inFileName, [sectors[1]], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeRhoP"] = fixedShapeRhoP
#	print "Finished with fixed shape rhoP"

#	print "Starting with fixed shape rhoF"
#	fixedShapeRhoF = doFixedShapes(inFileName, [sectors[2]], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeRhoF"] = fixedShapeRhoF
#	print "Finished with fixed shape rhoF"

#	print "Starting with fixed shape bothRho"
#	fixedShapeBothRho = doFixedShapes(inFileName, sectors[1:3], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeBothRho"] = fixedShapeBothRho
#	print "Finished with fixed shape bothRho"

#	print "Starting with fixed shape f2"
#	fixedShapeF2 = doFixedShapes(inFileName, [sectors[3]], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeF2"] = fixedShapeF2
#	print "Finished with fixed shape f2"

	print "Starting with fixed shapes"
	fixedShapes = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"
	fullSig = fixedShapes.getZeroModeSignature()

##	print "Starting with phase"
##	fitPiPiSshape = doF0phase(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave)
##	allMethods["pipiS"] = fitPiPiSshape
##	print "Finished with phase"

#	print "Starting with fitting rhoP"
#	fitRhoP = doFitRhoP(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fitRhoP"] = fitRhoP
#	print "Finished with fitting rhoP"

#	print "Starting with fitting rhoF"
#	fitRhoF = doFitRhoF(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fitRhoF"] = fitRhoF
#	print "Finished with fitting rhoF"

	print "Starting with fitting bothRho"
	fitBothRho = doFitBothRho(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fitBothRho"] = fitBothRho
	print "Finished with fitting bothRho"

	print "Starting with fitting f2"
	fitF2 = doFitF2(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fitF2"] = fitF2
	print "Finished with fitting f2"

	if stopBin - startBin > 1:
		print "Starting with smooth"
		smooth = doSmooth(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
		allMethods["smooth"] = smooth
		print "Finished with smooth"


	if "fixedShapeRhoF" in allMethods:
		allMethods["fixedShapeRhoF"].setZeroModeSignature(fullSig,1)
	if "fitRhoF" in allMethods:
		allMethods["fitRhoF"].setZeroModeSignature(fullSig,1)

	diffsFull, resolvedWA, nonResolvedWA, comps, resDiffs, nonResDiffs, resolvedDiffsFull,noCorrDiffs = cu.doAllComparisons(allMethods, startBin, methodBinRanges)
#	print resolvedDiffsFull
	from math import isnan
	for pair in resolvedDiffsFull:
		with  modernplotting.toolkit.PdfWriter('./resolvedDiffPlots/2mp_'+pair[0]+"_"+pair[1]+".pdf") as pdfOutput:
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

	with  modernplotting.toolkit.PdfWriter("compositions_2mp_MC.pdf") as pdfOutput:
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

	style.titleRight = r"$2^{-+}0^+$"
	style.titleLeft  = r"Monte Carlo"

	with modernplotting.toolkit.PdfWriter("studies_2mp.pdf") as pdfOutput:
		plot = style.getPlot2D()
		plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(len(studyList)+2)])
		plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(len(studyList))])
		studyPlotter.makeValuePlot(plot, hist)
		
		plot.axes.set_yticklabels(axolotl)
		axolotl.append(unCorrected_string)
		axolotl.append(weightedAVG_string)
		plot.axes.set_xticklabels(axolotl, rotation = 90)
		plot.setZlim((0.,1.))

		pdfOutput.savefigAndClose()

	with open("studies_2mp.txt",'w') as out:
		for axl in axolotl:
			out.write(axl + ' ')
		out.write("\n")
		for i in range(hist.GetNbinsX()):
			for j in range(hist.GetNbinsY()):
				out.write(str(hist.GetBinContent(i+1, j+1)) + ' ')
			out.write('\n')

	return

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
				fileName = "./collectedMethods/"+stu+"_"+sect+"_2mpMC_"+str(bin)
				if not (sect,bin) in fileNames:
					fileNames[sect,bin] = []
				fileNames[sect,bin].append(fileName)
				rv.writeAmplFiles(bin, fileName = fileName)

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
		allMethods['fixedShapes'].removeZeroModeFromComa()
		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = allMethods['fixedShapes'].produceResultViewer(resolvedWA,s, noRun = True, plotTheory = True)
		rv.writeBinToPdf(startBin, stdCmd = ["./comparisonResultsData/" + sect + "_MC_2D_"+str(tBin)+".pdf", "", [], "", []])
		for b in range(startBin, stopBin):
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			rv.writeBinToPdf(b, stdCmd = ["", "./comparisonResults/" + sect + "_MC_intens_"+str(b)+".pdf", intensNames,  "./comparisonResults/" + sect + "_MC_argand_"+str(b)+".pdf", argandNames])
	print studyList


if __name__ == "__main__":
	main()
