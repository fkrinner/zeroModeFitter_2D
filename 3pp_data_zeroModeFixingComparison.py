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
from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0,mRho3,Grho3
from rootfabi import root_open

import consistencyUtils as cu
import LaTeX_strings

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from LaTeX_strings import unCorrected_string, weightedAVG_string

def doFitRho(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {},referenceWave = "", writeResultToFile = None):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 2, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"3++0+[pi,pi]1--PiD":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
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
			resultString = str(tBins[0])+ " 666. " + str(rhoMass.value) + ' ' + str(rhoMass.error) + ' ' + str(rhoWidth.value) + ' ' + str(rhoWidth.error) + "\n"
			outFile.write(resultString)
	return fitRho

def doFitF2(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {},referenceWave = "", writeResultToFile = None):
	f2Mass  = ptc.parameter( mRho, "f2Mass" )
	f2Width = ptc.parameter( Grho , "f2Width")
	f2 = ptc.relativisticBreitWigner([f2Mass,f2Width], mPi, mPi, mPi, 1, 2, False)
	fitF2 = amplitudeAnalysis(inFileName, sectors, {"3++0+[pi,pi]2++PiP":[f2]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
	fitF2.loadData(referenceWave = referenceWave)
	fitF2.finishModelSetup()
	fitF2.fitShapeParameters()
	fitF2.calculateNonShapeParameters()
	fitF2.mode = AMPL
#	fitF2.removeGlobalPhaseFromComa()
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(f2Mass.value) + ' ' + str(f2Mass.error) + ' ' + str(f2Width.value) + ' ' + str(f2Width.error) + "\n"
			outFile.write(resultString)
	return fitF2

def doFitRho3(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {},referenceWave = "", writeResultToFile = None):
	rho3Mass  = ptc.parameter( mRho3, "rhoMass" )
	rho3Width = ptc.parameter( Grho3, "rhoWidth")
	rho3 = ptc.relativisticBreitWigner([rho3Mass,rho3Width], mPi, mPi, mPi, 3, 0, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"3++0+[pi,pi]3--PiS":[rho3]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, zeroFileName = zeroFileName)
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
			resultString = str(tBins[0])+ " 666. " + str(rho3Mass.value) + ' ' + str(rho3Mass.error) + ' ' + str(rho3Width.value) + ' ' + str(rho3Width.error) + "\n"
			outFile.write(resultString)
	return fitRho


alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def doFixedShapes(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {},referenceWave = ""):
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

def doSmooth(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {},referenceWave = ""):
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
	inFileName   = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_3pp.root"
	zeroFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_3pp.root"

	sectors          = ["3++0+[pi,pi]1--PiD", "3++0+[pi,pi]3--PiS"]
	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	tBins          = [tBin]
	sectorRangeMap = {}
	rhoRangeAdder  = ''
	startBin       = 11
	stopBin        = 50
	restrictRho    = False
	if restrictRho:
		rhoRange       = 1.2
		sectorRangeMap = {"3++0+[pi,pi]1--PiD":(0.,rhoRange)}
		rhoRangeAdder  = "_range"+str(rhoRange)

	methodBinRanges = {
	                   "fitF2"        : (22, 50),
	                   "fitRho3"      : (32, 50),
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
	                      "fitRho"            : r"$\text{fit}_\rho$",
	                      "fitRhoP"           : r"$\text{fit}_\rho^{P}$",
	                      "fitRhoF"           : r"$\text{fit}_\rho^{F}$",
	                      "fitBothRho"        : r"$\text{fit}_\rho^{2}$",
	                      "fitBothF2"         : r"$\text{fit}_{f_2}^{2}$",
	                      "fitF2"             : r"$\text{fit}_{f_2}$",
	                      "fitRho3"           : r"$\text{fit}_{\rho_3}$",
	                      "smooth"            : r"smooth"}



	print "Starting with fixed shapes"
	fixedShapes = doFixedShapes(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, sectorRangeMap = sectorRangeMap)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"
	fullSig = fixedShapes.getZeroModeSignature()


	print "Starting with fitting rho"
	fitRho = doFitRho(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "./rhoMassesAndWidths_3pp_global"+rhoRangeAdder+".dat", sectorRangeMap = sectorRangeMap)
	allMethods["fitRho"] = fitRho
	print "Finished with fitting rho"

	print "Starting with fitting rho3"
	fitRho3 = doFitRho3(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "./rho3MassesAndWidths_3pp_global.dat", sectorRangeMap = sectorRangeMap)
	allMethods["fitRho3"] = fitRho3
	print "Finished with fitting rho3"

	if stopBin - startBin > 1:
		print "Starting with smooth"
		smooth = doSmooth(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave, sectorRangeMap = sectorRangeMap)
		allMethods["smooth"] = smooth
		print "Finished with smooth"

	if "fixedShapeF0" in allMethods:
		allMethods["fixedShapeF0"].setZeroModeSignature(fullSig,1)

	diffsFull, resolvedWA, nonResolvedWA, comps, resDiffs, nonResDiffs, resolvedDiffsFull,noCorrDiffs = cu.doAllComparisons(allMethods, startBin, methodBinRanges)
#	print resolvedDiffsFull
	from math import isnan
	for pair in resolvedDiffsFull:
		with  modernplotting.toolkit.PdfWriter('./resolvedDiffPlots/bigger3pp_'+pair[0]+"_"+pair[1]+"_"+str(tBin)+".pdf") as pdfOutput:
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

	style.titleRight = r"$3^{++}0^+$"
	style.titleLeft  = LaTeX_strings.tBins[tBin]

#	with  modernplotting.toolkit.PdfWriter("compositions_bigger3pp_"+str(tBin)+".pdf") as pdfOutput:
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
#		pdfOutput.savefigAndClose()



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

	with modernplotting.toolkit.PdfWriter("studies_bigger3pp_"+str(tBin)+".pdf") as pdfOutput:
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

	with open("studies_bigger3pp_"+str(tBin)+".txt",'w') as out:
		for axl in axolotl:
			out.write(axl + ' ')
		out.write("\n")
		for i in range(hist.GetNbinsX()):
			for j in range(hist.GetNbinsY()):
				out.write(str(hist.GetBinContent(i+1, j+1)) + ' ')
			out.write('\n')

	doRhoFits  = False
	doF2Fits   = False
	doRho3Fits = True
	binWiseFit = True

	if doRhoFits:
		with open("rhoMassesAndWidths_3pp_"+str(tBin)+rhoRangeAdder+".dat",'w') as out:
			if binWiseFit:
				bins = range(stopBin-startBin)
			else:
				bins = [1]
			for i in bins:
				binIndex = i+startBin
				out.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				while True:
					try:
						if binWiseFit:
							fitRange = [i]
						else:
							fitRange = range(stopBin-startBin)

						x,err,c2,ndf = fitRho.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0], fitRange, zeroModeParameters = resolvedWA)
						break
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							print "Too many failed attempts in bin "+str(i)+": "+str(exceptCount)
#							raise Exception
							x, err = [0.,0.],[0.,0.]
							c2, ndf = 0.,1
				out.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				out.write(' '+ str(c2/ndf) +'\n')

			with open("3pp_rho_cpls_"+str(tBin)+".dat",'w') as outFileCpl:
				fitRho.calculateNonShapeParameters()
				cpl, hess = fitRho.getCouplingParameters()
				outFileCpl.write(str(cpl))
				outFileCpl.write("\n")
				outFileCpl.write(str(hess))

	if doF2Fits:
		fitF2 = doFitF2(inFileName,zeroFileName, ["3++0+[pi,pi]2++PiP"], startBin, stopBin, tBins)
		with open("f2MassesAndWidths_3pp_"+str(tBin)+".dat",'w') as out:
			if binWiseFit:
				bins = range(stopBin-startBin)
			else:
				bins = [1]
			for i in bins:
				binIndex = i+startBin
				out.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				while True:
					try:
						if binWiseFit:
							fitRange = [i]
						else:
							fitRange = range(stopBin-startBin)
						x,err,c2,ndf = fitF2.fitShapeParametersForBinRange([mF2+startValueOffset,GF2+startValueOffset], [0],fitRange, zeroModeParameters = [[[] for _ in range(stopBin-startBin)]])
						break
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							print "Too many failed attempts in bin "+str(i)+": "+str(exceptCount)
#							raise Exception
							x, err = [0.,0.],[0.,0.]
							c2, ndf = 0.,1
							break


				out.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				out.write(' '+ str(c2/ndf) +'\n')
	if doRho3Fits:
		with open("rho3MassesAndWidths_3pp_"+str(tBin)+".dat",'w') as out:
			if binWiseFit:
				bins = range(stopBin-startBin)
			else:
				bins = [1]
			for i in bins:
				binIndex = i+startBin
				out.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				while True:
					try:
						if binWiseFit:
							fitRange = [i]
						else:
							fitRange = range(stopBin-startBin)

						x,err,c2,ndf = fitRho3.fitShapeParametersForBinRange([mRho3+startValueOffset,Grho3+startValueOffset], [0],fitRange, zeroModeParameters = resolvedWA)
						break
					except:
						print "Fitter exception encountered"
						startValueOffset += 0.001
						exceptCount      += 1	
						if exceptCount > 3:
							print "Too many failed attempts in bin "+str(i)+": "+str(exceptCount)
#							raise Exception
							x, err = [0.,0.],[0.,0.]
							c2, ndf = 0.,1
							break

				fitRho3.calculateNonShapeParametersForZeroModeParameters(resolvedWA)
				rhoRV = fitRho3.produceResultViewer(resolvedWA,"3++0+[pi,pi]3--PiS", noRun = True, plotTheory = True)
				rhoRV.plotData = True
				plotNameBase = "./rho3FitPlots/3pp0p3mmPiS_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
				rhoRV.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])

				with open("3pp_rho3_cpls_"+str(tBin)+".dat",'w') as outFileCpl:
					fitRho3.calculateNonShapeParameters()
					cpl, hess = fitRho3.getCouplingParameters()
					outFileCpl.write(str(cpl))
					outFileCpl.write("\n")
					outFileCpl.write(str(hess))
				out.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				out.write(' '+ str(c2/ndf) +'\n')
	if doRhoFits or doF2Fits or doRho3Fits:			
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
				fileName = "./collectedMethods/"+stu+"_"+sect+"_bigger3pp_"+str(bin)
				if not (sect,bin) in fileNames:
					fileNames[sect,bin] = []
				fileNames[sect,bin].append(fileName)
				rv.writeAmplFiles(bin, fileName = fileName)

#	totalHists = fixedShapes.getTotalHists(resolvedWA)
#	with root_open("./totals_3pp.root", "UPDATE") as out:
#		for t in totalHists:
#			for m in t:
#				m.Write()

	plotFolder      = "./comparisonResults3pp/"
	plotFolderScale = "./comparisonResults3pp_scale/"

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
		allMethods['fixedShapes'].removeZeroModeFromComa()
		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rvNC = allMethods['fixedShapes'].produceResultViewer(cloneZeros(resolvedWA),s, noRun = True, plotTheory = True)
		rvNC.writeBinToPdf(startBin, stdCmd = [plotFolder + sect + "_data_2D_"+str(tBin)+"_noCorr.pdf", "", [], "", []])
		rv = allMethods['fixedShapes'].produceResultViewer(resolvedWA,s, noRun = True, plotTheory = True)
		rvSC = allMethods['fixedShapes'].produceResultViewer(resolvedWA,s, noRun = True, plotTheory = True)
		rv.writeBinToPdf(startBin, stdCmd = [plotFolder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
		rvSC.scaleTo = "maxCorrData"
		for b in range(startBin, stopBin):
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			rv.writeBinToPdf(b, stdCmd = ["", plotFolder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  plotFolder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
			rvSC.writeBinToPdf(b, stdCmd = ["", plotFolderScale + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  plotFolderScale + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
	print studyList


if __name__ == "__main__":
	main()
