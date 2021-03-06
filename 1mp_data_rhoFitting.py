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
import numpy.linalg as la
from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0
from rootfabi import root_open
import LaTeX_strings

import consistencyUtils as cu

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from LaTeX_strings import unCorrected_string, weightedAVG_string
def doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None, acv = None):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"1-+1+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(loadIntegrals = True, referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	if writeResultToFile is not None:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(rhoMass.value) + ' ' + str(rhoMass.error) + ' ' + str(rhoWidth.value) + ' ' + str(rhoWidth.error) + "\n"
			outFile.write(resultString)
	return fitRho

def doFitRhoRadius(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None, acv = None):
	rhoMass   = ptc.parameter( mRho, "rhoMass" )
	rhoWidth  = ptc.parameter( Grho, "rhoWidth")
	rhoRadius = ptc.parameter(  Pr0, "rRho"    )
	pi1radius = ptc.parameter(  100., "rPi1"    )
	
	rhoMass.lock   = True
	rhoWidth.lock  = True
#	pi1radius.lock = True
	rhoRadius.lock = True

	rho = ptc.relativisticBreitWigner([rhoMass, rhoWidth, pi1radius, rhoRadius], mPi, mPi, mPi, 1, 1, fitPr = True)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"1-+1+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData( loadIntegrals = True, referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	if writeResultToFile is not None:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0]) + " 666. " + str(rhoRadius.value) + "\n"
#	#	#	resultString = str(tBins[0])+ " 666. "  + str(rhoMass.value) + ' ' + str(rhoMass.error) + ' ' + str(rhoWidth.value) + ' ' + str(rhoWidth.error) + "\n"
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

	inFileName    = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_exotic.root"
	sectors          = ["1-+1+[pi,pi]1--PiP"]
	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	tBins            = [tBin]

	startBin         = 13
	stopBin          = 50

	acv = 130. # artificial coma value
#	acv = 1.3e-5 # artificial coma value


#	startBin = 26
#	stopBin  = 28

	methodBinRanges = {
	                   "fitF2"        : (22, 50),
	                   "fitRhoF"      : (22, 50),
	                   "fixedShapeF2" : (22, 50)}
#	methodBinRanges = {} # Override here 


#	sectorRangeMap = {"1-+1+[pi,pi]1--PiP":(0.,1.2)}
	sectorRangeMap = {"1-+1+[pi,pi]1--PiP":(0.,1.12)}	

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
	fixedShapes = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, acv = acv, sectorRangeMap = sectorRangeMap)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"
	fullSig = fixedShapes.getZeroModeSignature()


#	print "Starting with fitting rho"
#	fitRho = doFitRho(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "rhoMassesAndWidths_1-+1+1--_global.dat", sectorRangeMap = sectorRangeMap, acv = acv)
#	allMethods["fitRho"] = fitRho
#	print "Finished with fitting rho"

#	if stopBin - startBin > 1:
#		print "Starting with smooth"
#		smooth = doSmooth(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
#		allMethods["smooth"] = smooth
#		print "Finished with smooth"

	if "fixedShapeRhoF" in allMethods:
		allMethods["fixedShapeRhoF"].setZeroModeSignature(fullSig,1)
	if "fitRhoF" in allMethods:
		allMethods["fitRhoF"].setZeroModeSignature(fullSig,1)

#	diffsFull, resolvedWA, nonResolvedWA, comps, resDiffs, nonResDiffs, resolvedDiffsFull,noCorrDiffs = cu.doAllComparisons(allMethods, startBin, methodBinRanges, noBelowZero = True)

#	print resolvedDiffsFull
#	from math import isnan
#	for pair in resolvedDiffsFull:
#		with  modernplotting.toolkit.PdfWriter('./resolvedDiffPlots/1mp_'+pair[0]+"_"+pair[1]+"_"+str(tBin)+".pdf") as pdfOutput:
#			plot  = style.getPlot1D()
#			line  = [0.000000001]*len(resolvedDiffsFull[pair][0])
#			line2 = [0.000000001]*len(resolvedDiffsFull[pair][0])
#			one   = [1.]*len(resolvedDiffsFull[pair][0])
#			xAxis = [ .5 + 0.04*(startBin + i) for i in range(len(resolvedDiffsFull[pair][0]))]
#			for i,v in enumerate(resolvedDiffsFull[pair][0]):
#				if not isnan(v) and not v <= 0.:
#					line[i] = v
#				else:
#					line[i] = 0.000000001
#
#			if not pair[1] == "WAres" and not pair[1] == "WAnon":
#				for i,v in enumerate(resolvedDiffsFull[pair[1],pair[0]][0]):
#					if not isnan(v) and not v <= 0.:
#						line2[i] = v
#					else:
#						line2[i] = 0.000000001
#
#			plot.setYlog()
#			plot.plot(xAxis, line)
#			plot.plot(xAxis, one)
#			plot.plot(xAxis, line2)
#			plot.setYlim(0.00001, 10000)
#			plot.finishAndSaveAndClose(pdfOutput)

#	studyList = []
#	for m in allMethods:
#		studyList.append(m)
#	studyList.sort()

	style.titleRight = r"$1^{-+}1^+$"
	style.titleLeft  = LaTeX_strings.tBins[tBin]

#	with  modernplotting.toolkit.PdfWriter("compositions_1mp_"+str(tBin)+".pdf") as pdfOutput:
#		plot = style.getPlot1D()
#		for m in studyList:
#			line  = [0.]*len(comps[m][0])
#			xAxis = [ .5 + 0.04*(startBin + i) for i in range(len(comps[m][0]))]
#			break
#		count = 0
#		for m in studyList:
#			newLine = line[:]
#			for i in range(len(comps[m][0])):
#				newLine[i] += comps[m][0][i]
#			plot.axes.fill_between(xAxis, line, newLine, facecolor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, 0.1*count))
#			count += 1
#			line = newLine
#		plot.setYlim(0.,1.)
#		plot.setXlim(xAxis[0], xAxis[-1])
#		plot.finishAndSaveAndClose(pdfOutput)

#	hist = pyRootPwa.ROOT.TH2D("hist","hist", len(studyList)+2, 0, len(studyList)+2, len(studyList), 0, len(studyList))
#
#	for i,m in enumerate(studyList):
#		for j,n in enumerate(studyList):
#			hist.SetBinContent(i+1, j+1, diffsFull[n,m])
#	for i,m in enumerate(studyList):
#		hist.SetBinContent(len(studyList)+1, i+1, noCorrDiffs[m])
#		hist.SetBinContent(len(studyList)+2, i+1, resDiffs[m])
#
#	axolotl = []
#	for i,study in enumerate(studyList):
#		axolotl.append(shortlabels[study])
#		axolotl.append(alphabet[i])

#	with modernplotting.toolkit.PdfWriter("studies_1mp_"+str(tBin)+".pdf") as pdfOutput:
#		plot = style.getPlot2D()
#		plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(len(studyList)+2)])
#		plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(len(studyList)  )])
#		studyPlotter.makeValuePlot(plot, hist)
#		
#		plot.axes.set_yticklabels(axolotl)
#		axolotl.append(unCorrected_string)
#		axolotl.append(weightedAVG_string)
#		plot.axes.set_xticklabels(axolotl, rotation = 90)
#		plot.setZlim((0.,1.))
#
#		plot.finishAndSaveAndClose(pdfOutput)

#	with open("studies_1mp_"+str(tBin)+".txt",'w') as out:
#		for axl in axolotl:
#			out.write(axl + ' ')
#		out.write("\n")
#		for i in range(hist.GetNbinsX()):
#			for j in range(hist.GetNbinsY()):
#				out.write(str(hist.GetBinContent(i+1, j+1)) + ' ')
#			out.write('\n')

	doRhoFits = True
	writeCpls = False
	if writeCpls:
		outFileCpl = open("1mp_rho_cpls_"+str(tBin)+".dat",'w') 
	resolvedWA = allMethods['fixedShapes'].getZeroModeParametersForMode()
	doActuallyNotFit = False
	if doRhoFits:
		with open("pi1Radii_exotic_"+str(tBin)+".dat",'w') as out:
			rhoRadius = doFitRhoRadius(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = referenceWave, writeResultToFile = "rhoRadiusExotic_global.dat", acv = acv)
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				out.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				if doActuallyNotFit:
					print "The fit has been turned off, due to a workaround... :("
				else:
					exceptCount      = 0
					try:
						x,err,c2,ndf = rhoRadius.fitShapeParametersForBinRange([Pr0+ startValueOffset], [0],[i], zeroModeParameters = resolvedWA)
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							raise Exception("Too many failed attempts: "+str(exceptCount))
	
					out.write(str(x[0]) + ' ' + str(err[0]))
					out.write(' '+str(c2/ndf)+'\n')

				if writeCpls:
					rhoRadius.calculateNonShapeParametersForZeroModeParameters(resolvedWA)
					cpl, hess = rhoRadius.getCouplingParameters()
					hessInv = la.inv(hess[0][i])
					if not len(cpl[0][1]) == 2:
						raise IndexError("Parameter count not 2, change implementation")
					integral = rhoRadius.model[0][i].getIntegralForFunction(0, rhoRadius.model[0][i].funcs[0][0])
					outFileCpl.write(str(0.52 + binIndex*.04) + ' ' + str(cpl[0][i][0]**2 + cpl[0][i][0]**2) + ' ' + str(integral) + ' ')
					outFileCpl.write(str(cpl[0][i][0])     + ' ' + str(cpl[0][i][1])     + ' ')
					outFileCpl.write(str(hessInv[0][0]/2) + ' ' + str(hessInv[1][1]/2) + ' ' + str(hessInv[0][1]/2))
					outFileCpl.write("\n")

				makeRhoFitPlots = False
				if makeRhoFitPlots:
					rhoRadius.calculateNonShapeParametersForZeroModeParameters(resolvedWA)
					rhoRV = rhoRadius.produceResultViewer(resolvedWA,"1-+1+[pi,pi]1--PiP", noRun = True, plotTheory = True)
					rhoRV.plotData = True
					plotNameBase = "./rhoFitPlots/1mp1p1mmPiP_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
					rhoRV.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
		if writeCpls:
			outFileCpl.close()

		return
##### Writing starts here

#	fileNames = {}

#	fixedShapes.removeZeroModeFromComa()
##	binRange = {"1-+1+[pi,pi]1--PiP": (10,22)}
#	totalHists = fixedShapes.getTotalHists(resolvedWA, binRange = binRange)
#	with root_open("./totals_1mp_rhoRange.root", "UPDATE") as out:
#		for t in totalHists:
#			for m in t:
#				m.Write()
#	return

#	for stu in allMethods:
#		print "Writing for '" + stu + "'"
#		for s, sect in enumerate(allMethods[stu].sectors):
#			if stu == "pipiS":
#				rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, plotTheory = True)
#				rv.run()
#			rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, noRun = True)
#			for bin in range(startBin, stopBin):
#				fileName = "./collectedMethods/"+stu+"_"+sect+"_1mpData_"+str(bin)
#				if not (sect,bin) in fileNames:
#					fileNames[sect,bin] = []
#				fileNames[sect,bin].append(fileName)
#				rv.writeAmplFiles(bin, fileName = fileName)
	scalle = False
	if scalle:
		folder = "./comparisonResultsData_1mp_scale/"
	else:
		folder = "./comparisonResultsData_1mp/"

	fixedShapesFull = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, acv = acv)

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
#		allMethods['fixedShapes'].removeZeroModeFromComa()
#		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = fixedShapesFull.produceResultViewer(allMethods['fixedShapes'].getZeroModeParametersForMode(),s, noRun = True, plotTheory = True)
#		rv.plotData = False
		rv.writeBinToPdf(startBin, stdCmd = [ folder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
#		rv.labelPoints     = [0,10,15,20]
#		rv.makeLegend      = True
		if scalle:
			rv.scaleTo         = "maxCorrData"
#			rv.yAxisShift      = 300.
#			rv.tStringYpos     = 0.8
#			rv.topMarginIntens = 1.4
#			fakkkk             = 1.
		else:
			rv.plotData        = False
#			fakkkk             = .7	
#			rv.tStringYpos     = 0.865
#			rv.topMarginIntens = 1.3
#			rv.yAxisShift      = 100.
#		rv.shiftMap        = {0:(fakkkk*50.,fakkkk*-280.),10:(fakkkk*-420.,fakkkk*-50.), 15:(fakkkk*-420.,fakkkk*-30.), 20:(fakkkk*-50.,fakkkk*70.)}
		for b in range(startBin, stopBin):
#			if not bin == 27:
#				continue
#			intensNames = [name+".intens" for name in fileNames[sect,b]]
#			argandNames = [name+".argand" for name in fileNames[sect,b]]
#			intensNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/evalDima/dima.intens"]
#			argandNames = ["/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/evalDima/dima.argand"]
#			if not scalle:
			if True:
				intensNames = []
				argandNames = []
#			rv.plotData = True
			rv.writeBinToPdf(b, stdCmd = ["", folder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  folder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
			rv.wiriteReImToPdf(b, folder +  sect + "_data_<ri>_"+str(b)+"_"+str(tBin)+".pdf")


if __name__ == "__main__":
	main()
