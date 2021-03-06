import matplotlib
matplotlib.use("Agg")

from analysisClass import amplitudeAnalysis
import parameterTrackingParameterizations as ptc
import parameterizationClasses as pc
from fixedparameterizationPaths import getFileNameForSector
from modes import PHASE, AMPL, SMOOTH, NONE
from utils import sumUp, weightedSum, cloneZeros, checkLaTeX
import sys,os
import pyRootPwa
import numpy as np
from globalDefinitions import referenceWave, mPi, mK, mRho, Grho, mRhoPrime, GrhoPrime, mF0, g1, g2, mF2, GF2, Pr0, m1500, G1500
from rootfabi import root_open
import LaTeX_strings

import consistencyUtils as cu
from studyFileNames import fileNameMap

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from LaTeX_strings import unCorrected_string, weightedAVG_string
def doF0Fit(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", deg = 0, writeResultToFile = None):
	fixedAMPD_sub = pc.fixedParameterization("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_0mp0pSigmaPiS")

	f0Mass   = ptc.parameter(mF0,"f0Mass")
#	f0Mass.lock = True
	f0g1     = ptc.parameter(g1, "g1")
	f0g2     = ptc.parameter(g2, "g2")
#	f0g1.lock = True
#	f0g2.lock = True
	f0       = ptc.flatte([f0Mass, f0g1, f0g2], mPi, mPi, mPi, mK, mK, 0, 0, False)
	polyPars = []
	for i in range(deg):
		polyPars.append(ptc.parameter(0.1, "reC"+str(i+1)))
		polyPars.append(ptc.parameter(0.1, "imC"+str(i+1)))
	poly  = ptc.complexPolynomial(deg, polyPars)
	fitF0 = amplitudeAnalysis(inFileName, sectors, {"0-+0+[pi,pi]0++PiS":[f0, fixedAMPD_sub]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitF0.loadData(referenceWave = referenceWave)
	fitF0.finishModelSetup()
	fitF0.fitShapeParameters()
	fitF0.calculateNonShapeParameters()
	fitF0.mode = AMPL
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(f0Mass.value) + ' ' + str(f0Mass.error) + ' ' + str(f0g1.value) + ' ' + str(f0g1.error) + ' ' + str(f0g2.value) + ' ' + str(f0g2.error) 
			for i in range(deg):
				rasultString += ' ' + str(polyPars[i].value) + ' ' + str(polyPars[i].error)
			resultString +=  "\n"
			outFile.write(resultString)
	return fitF0

def doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}, referenceWave = "", writeResultToFile = None):
	rhoMass  = ptc.parameter( mRho,  "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"0-+0+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
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

def doFit1500(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]0++PiS": (1.15, 1.75)}, referenceWave = "", writeResultToFile = None):
	fixedAMPD_sub = pc.fixedParameterization("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_0mp0pSigmaPiS")
	const = pc.constant()
	linea = pc.linear()
	f0mass  = ptc.parameter( m1500,  "1500mass" )
	f0width = ptc.parameter( G1500 , "1500width")
	f0_1500 = ptc.relativisticBreitWigner([f0mass,f0width], mPi, mPi, mPi, 0, 0, False)
	fit1500 = amplitudeAnalysis(inFileName, sectors, {"0-+0+[pi,pi]0++PiS":[f0_1500, fixedAMPD_sub]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
#	fit1500 = amplitudeAnalysis(inFileName, sectors, {"0-+0+[pi,pi]0++PiS":[f0_1500, const, linea]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fit1500.loadData(referenceWave = referenceWave)
	fit1500.finishModelSetup()
	fit1500.fitShapeParameters()
	fit1500.calculateNonShapeParameters()
	fit1500.mode = AMPL
#	fit1500.removeGlobalPhaseFromComa()
	if writeResultToFile:
		with open(writeResultToFile, 'a') as outFile:
			if len(tBins) > 1:
				raise ValueError("More than one t' bin not supported")
			resultString = str(tBins[0])+ " 666. " + str(f0mass.value) + ' ' + str(f0mass.error) + ' ' + str(f0width.value) + ' ' + str(f0width.error) + "\n"
			outFile.write(resultString)
	return fit1500

def doF0phase(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]0++PiS" : (0.34, 2*mK - .04)}, referenceWave = ""):
	pipiSfileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_0mp0pSigmaPiS"
	pipiSw = pc.fixedParameterization(pipiSfileName, polynomialDegree  = 0, complexPolynomial = False)
	waveModel = {"0-+0+[pi,pi]0++PiS": [pipiSw]}
	fitPiPiSshape = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitPiPiSshape.loadData(referenceWave = referenceWave)
	fitPiPiSshape.finishModelSetup()
	fitPiPiSshape.phaseFit()
	fitPiPiSshape.mode = PHASE
#	fitPiPiSshape.removeGlobalPhaseFromComa()
	fitPiPiSshape.unifyComa()
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
	fixedShapes.loadData(loadIntegrals = True,  referenceWave = referenceWave)
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

	sectors          = ["0-+0+[pi,pi]0++PiS", "0-+0+[pi,pi]1--PiP"]
#	rhoRange         = 1.2
	rhoRange         = None
	if not rhoRange:
		sectorRangeMap = {}
		rhoRangeString = ""
	else:
		sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(0.,rhoRange)}
		rhoRangeString = "_range"+str(rhoRange)


	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))
	if len(sys.argv) > 2:
		study      = sys.argv[2]
		studyAdder = "_"+study
	else:
		study      = "std11"
		studyAdder = ""
	print "Study: "+study

	referenceWave    = ""

	inFileName       = fileNameMap[study]
	tBins            = [tBin]

	startBin         = 11
	stopBin          = 50

#	startBin         = 32
#	stopBin          = 33

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

	fit1500 = doFit1500(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave, writeResultToFile = "0mp_f0_1500_massesAndWidths_global.dat")

	print "Start with fixed shape f0"
	fixedShapeF0 = doFixedShapes(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"

#	print "Start with fixed shape rho"
#	fixedShapeRho = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fixedShapeRho"] = fixedShapeRho
#	print "Finished with fixed shape rho"

#	print "Start with restricted rho (1 Gamma)"
#	fixedShapeRho1G = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho - Grho, mRho+Grho)}, referenceWave = referenceWave)
#	allMethods["fixedShapeRho1G"] = fixedShapeRho1G
#	print "Finished with restricted rho (1 Gamma)"

#	print "Start with restricted rho (2 Gammas)"
#	fixedShapeRho2G = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho -2*Grho, mRho+2*Grho)}, referenceWave = referenceWave)
#	allMethods["fixedShapeRho2G"] = fixedShapeRho2G
#	print "Finished with restricted rho (2 Gammas)"

	print "Start with fixed shapes"
	fixedShapes = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"


#	print "Start with phase"
#	fitPiPiSshape = doF0phase(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["pipiS"] = fitPiPiSshape
#	print "Finished with phase"

	print "Start with fitting rho"
	fitRho = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap, referenceWave = referenceWave, writeResultToFile = "rhoMassesAndWidths_0-+0+1--_global.dat")
#	allMethods["fitRho"] = fitRho
	print "Finished with fitting rho"

#	print "Start with fitting restricted rho (1 Gamma)"
#	fitRho1G = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho - Grho, mRho+Grho)}, referenceWave = referenceWave)
#	allMethods["fitRho1G"] = fitRho1G
#	print "Finished with fitting restricted rho (1 Gamma)"

#	print "Start with fitting restricted rho (2 Gammas)"
#	fitRho2G = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho -2*Grho, mRho+2*Grho)}, referenceWave = referenceWave)
#	allMethods["fitRho2G"] = fitRho2G
#	print "Finished with fitting restricted rho (2 Gammas)"

	if stopBin - startBin > 1:
		print "Start with smooth"
		smooth = doSmooth(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
		allMethods["smooth"] = smooth
		print "Finished with smooth"


	ndfs   = {}
	params = {}
	for m in allMethods:
		ndfs[m]=  allMethods[m].getNDFforMode()
		params[m] = allMethods[m].getZeroModeParametersForMode()
		print m,sumUp(allMethods[m].evaluateZeroModeParametersForMode(params[m])).real/ndfs[m]
	diffs = cu.getmBinResolvedDiffs(allMethods)
	comps = cu.getCompositions(diffs)
#	with  modernplotting.toolkit.PdfWriter("compositions_0mp_data"+str(tBin)+studyAdder+".pdf") as pdfOutput:
#		plot = style.getPlot1D()
#		for m in comps:
#			line  = [0.]*len(comps[m][0])
#			xAxis = [ .5 + 0.04*(startBin + i) for i in range(len(comps[m][0]))]
#			break
#		count = 0
#		for m in comps:
#			print m
#
#			newLine = line[:]
#			for i in range(len(comps[m][0])):
#				newLine[i] += comps[m][0][i]
#			plot.axes.fill_between(xAxis, line, newLine, facecolor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, 0.1*count))
#			count += 1
#			line = newLine
#		plot.setYlim(0.,1.)
#		plot.setXlim(xAxis[0], xAxis[-1])
#		pdfOutput.savefigAndClose()
	studyList = []
	for m in allMethods:
		studyList.append(m)
	studyList.sort()

	selfEvals = {}
	for m in allMethods:
		selfEvals[m] = sumUp(allMethods[m].evaluateResolvedZeroModeParametersForMode(params[m])).real

	hist = pyRootPwa.ROOT.TH2D("hist","hist", len(params)+2, 0, len(params)+2, len(params), 0, len(params))

	cumulWeights = {}
#	resolvedWeightedSum = [[]] # Assumes one t' bin
#	for i in range(stopBin - startBin):
#		dim = len(params[m][0][i])
#		prrs = [0.] * dim
#		for m in params:
#			weight = comps[m][0][i]
#			if not m in cumulWeights:
#				cumulWeights[m] = 0.
#			cumulWeights[m] += weight
#			for j in range(dim):
##				print dim,i,j,params[m][0]
##				print prrs
#				if weight < 0.:
#					raise ValueError("Smaller than zero weight encountered...")
#				prrs[j] += weight * params[m][0][i][j]
#		resolvedWeightedSum[0].append(prrs)
#
#	evals = {}
#	for i,m in enumerate(studyList):
##		print "-------------------------------"
#		for j,n in enumerate(studyList):
#			evl = sumUp(allMethods[n].evaluateResolvedZeroModeParametersForMode(params[m])).real
#			evals[n,m] = evl
#			diff = (evl-selfEvals[n])/selfEvals[n]
#
##			allMethods["fixedShapes"].removeZeroModeFromComa()
##			print "------------------------------------IN---------------------------------"
##			print params[m], params[n]
##			diff = sumUp(allMethods["fixedShapes"].compareTwoZeroModeCorrections(params[m], params[n]))
##			print diff
##			print "------------------------------------OUT---------------------------------"
##			print m,'in',n,":",diff
#			hist.SetBinContent(i+1, j+1, diff)
##	return 
#	weightedSum = weightedParametersSum(evals, selfEvals, params)
#	for i,m in enumerate(studyList):
#		evl = sumUp(allMethods[m].evaluateZeroModeParametersForMode(cloneZeros(weightedSum))).real
#		diff = (evl - selfEvals[m])/selfEvals[m]
#		evl2 = sumUp(allMethods[m].evaluateZeroModeParametersForMode(resolvedWeightedSum)).real
#		diff2 = (evl2 - selfEvals[m])/selfEvals[m]
#
#		print m,diff,";:;:;;>>>??"
#		hist.SetBinContent(len(studyList)+1, i+1, diff)
#		hist.SetBinContent(len(studyList)+2, i+1, diff2)
#
#
#
#	axolotl = []
#	for i,study in enumerate(studyList):
#		axolotl.append(shortlabels[study])
#		axolotl.append(alphabet[i])

	style.titleRight = r"$0^{-+}0^+$"
	style.titleLeft  = LaTeX_strings.tBins[tBin]

#	with modernplotting.toolkit.PdfWriter("studies_0mp_data"+str(tBin)+studyAdder+".pdf") as pdfOutput:
#		plot = style.getPlot2D()
#		plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(len(studyList)+2)])
#		plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(len(studyList))])
#		studyPlotter.makeValuePlot(plot, hist)

#		plot.axes.set_yticklabels(axolotl)
#		axolotl.append(unCorrected_string)
#		axolotl.append(weightedAVG_string)
#		plot.axes.set_xticklabels(axolotl, rotation = 90)
#		plot.setZlim((0.,1.))

#		pdfOutput.savefigAndClose()


#	with open("studies_0mp_data"+str(tBin)+studyAdder+".txt",'w') as out:
#		for axl in axolotl:
#			out.write(axl + ' ')
#		out.write("\n")
#		for i in range(hist.GetNbinsX()):
#			for j in range(hist.GetNbinsY()):
#				out.write(str(hist.GetBinContent(i+1, j+1)) + ' ')
#			out.write('\n')

##### Writing starts here
	doF0fitGlobal = False
	doF0Fits      = False
	doFits1500    = False
	doRhoFits     = False
	if doF0fitGlobal:
		deg = 0
#		f0fit = doF0Fit(inFileName, sectors[:1], startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]0++PiS":(0.9,1.1)}, referenceWave = referenceWave, deg = deg)
		f0fit = doF0Fit(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave, deg = deg)
		startPars = [mF0,g1,g2]
		for _ in range(deg):
			startPars.append(0.)
			startPars.append(0.)
		x,err = f0fit.fitShapeParametersForBinRange(startPars,[0],range(stopBin-startBin), zeroModeParameters = resolvedWeightedSum)
	######	x,err = f0fit.fitShapeParametersForBinRange([mF0],[0],range(stopBin-startBin), zeroModeParameters = resolvedWeightedSum)
		print x
		f0fit.setShapeParameters(x,err,resolvedWeightedSum)
		s = 0 # only one sector??
		rv = f0fit.produceResultViewer(resolvedWeightedSum,s, noRun = True, plotTheory = True)
		for b in range(startBin, stopBin):
			intensName = "./f0fits/0mp_intens_"+str(b)+"_"+str(tBin)+".pdf"
			argandName = "./f0fits/0mp_argand_"+str(b)+"_"+str(tBin)+".pdf"
			rv.writeBinToPdf(b, stdCmd = ["", intensName, [],  argandName, []])

	if doF0Fits:
		f0fit = doF0Fit(inFileName, sectors[:1], startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]0++PiS":(0.9,1.1)}, referenceWave = referenceWave, writeResultToFile = "f0MassesAndWidths_0-+0+0++_global.dat")
		with open("./f0MassesAndWidths_0mp_"+str(tBin)+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
#				try:
				if True:
					startPars = [mF0+startValueOffset,g1+startValueOffset,g2+startValueOffset]
#					for _ in range(deg):
#						startPars.append(startValueOffset)
#						startPars.append(startValueOffset)
					x,err,c2,ndf = f0fit.fitShapeParametersForBinRange(startPars, [0],[i], zeroModeParameters = resolvedWeightedSum)
#				except:
#					print "Fitter exception encountered"
#					startValueOffset += 0.001
#					exceptCount      += 1
#					if exceptCount > 3:
#						raise Exception("Too many failed attempts: "+str(exceptCount))

				f0fit.calculateNonShapeParametersForZeroModeParameters(resolvedWeightedSum)
				RV980 = f0fit.produceResultViewer(resolvedWeightedSum,"0-+0+[pi,pi]0++PiS", noRun = True, plotTheory = True)
				RV980.plotData = True
				plotNameBase = "./980FitPlots/0mp0p0ppPiS_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
				RV980.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
				if binIndex == 32:
					RV980.writeToRootFile("980.root")
				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1])+ ' ' + str(x[2]) + ' ' + str(err[2]))
#				for i in range(deg):
#					outFile.write( ' ' + str(x[3+i]) + ' ' + str(err[3+i]) + ' ' + str(x[4+i]) + ' ' +str(err[4+i]))
				outFile.write(' ' +str(c2/ndf)+'\n')

	if doFits1500:
		with open("./0mp_f0_1500massesAndWidths_"+str(tBin)+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i + startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				if True:
					startPars = [m1500+startValueOffset,G1500+startValueOffset]
#					for _ in range(deg):
#						startPars.append(startValueOffset)
#						startPars.append(startValueOffset)
					x,err,c2,ndf = fit1500.fitShapeParametersForBinRange(startPars, [0],[i], zeroModeParameters = resolvedWeightedSum)
#				except:
#					print "Fitter exception encountered"
#					startValueOffset += 0.001
#					exceptCount      += 1
#					if exceptCount > 3:
#						raise Exception("Too many failed attempts: "+str(exceptCount))

				
				fit1500.calculateNonShapeParametersForZeroModeParameters(resolvedWeightedSum)
				RV1500 = fit1500.produceResultViewer(resolvedWeightedSum,"0-+0+[pi,pi]0++PiS", noRun = True, plotTheory = True)
				RV1500.plotData = True
				plotNameBase = "./1500FitPlots/0mp0p0ppPiS_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
				RV1500.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])
				if binIndex == 32:
					RV1500.writeToRootFile("1500.root")


				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
#				for i in range(deg):
#					outFile.write( ' ' + str(x[3+i]) + ' ' + str(err[3+i]) + ' ' + str(x[4+i]) + ' ' +str(err[4+i]))
				if not ndf == 0:
					outFile.write(' ' +str(c2/ndf)+'\n')
				else:
					outFile.write(' noNDF \n')

	if doRhoFits:
		cplFileName = "0mp_rho_cpls_"+str(tBin)+".dat"
		try:
			os.remove(cplFileName)
		except:
			pass
		with open("./rhoMassesAndWidths_0mp_"+str(tBin)+rhoRangeString+".dat",'w') as outFile:
			for i in range(stopBin-startBin):
				binIndex = i+startBin
				outFile.write(str(binIndex)+' '+str(0.52 + 0.04*binIndex)+' ')
				startValueOffset = 0.00
				exceptCount      = 0
				try:
					x,err,c2,ndf = fitRho.fitShapeParametersForBinRange([mRho+startValueOffset,Grho+startValueOffset], [0],[i], zeroModeParameters = resolvedWeightedSum)
				except:
					print "Fitter exception encountered"
					startValueOffset += 0.001
					exceptCount      += 1
					if exceptCount > 3:
						raise Exception("Too many failed attempts: "+str(exceptCount))

				outFile.write(str(x[0]) + ' ' + str(err[0]) + ' ' + str(x[1]) + ' ' + str(err[1]))
				outFile.write(' ' +str(c2/ndf)+'\n')

				fitRho.calculateNonShapeParametersForZeroModeParameters(resolvedWeightedSum)
				rhoRV = fitRho.produceResultViewer(resolvedWeightedSum,"0-+0+[pi,pi]1--PiP", noRun = True, plotTheory = True)
				rhoRV.plotData = True
				plotNameBase = "./rhoFitPlots/0mp0p1mmPiP_<mode>_"+str(binIndex)+"_"+str(tBin)+".pdf"
				rhoRV.writeBinToPdf(binIndex, stdCmd = ["", plotNameBase.replace("<mode>","intens"), [],  plotNameBase.replace("<mode>","argand"), []])

#				with open(cplFileName,'a') as outFileCpl:
#					fitRho.calculateNonShapeParameters()
#					cpl, hess = fitRho.getCouplingParameters()
#					outFileCpl.write(str(cpl[0][i]))
#					outFileCpl.write(", ")
#					outFileCpl.write(str(hess[0][i]))
#					outFileCpl.write("\n")

	if doF0fitGlobal or doF0Fits or doFits1500 or doRhoFits:
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


	resolvedWeightedSum = fixedShapes.getZeroModeParametersForMode()

	totalHists = fixedShapes.getTotalHists(resolvedWeightedSum)
#	studyAdder = "_noCorr"

	with root_open("./totals_0mp"+studyAdder+".root", "UPDATE") as out:
		for t in totalHists:
			for m in t:
				m.Write()
#	return
#	folder = "./comparisonResultsData"+studyAdder+"/"
	folder = "./"
	for s, sect in enumerate(allMethods['fixedShapes'].sectors):

#		allMethods['fixedShapes'].removeZeroModeFromComa()
#		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = allMethods['fixedShapes'].produceResultViewer(resolvedWeightedSum,s, noRun = True, plotTheory = True)
#		rv.titleFontSize = 11
		rv.showColorBar  = True
		rv.writeToRootFile("0mp_"+sect+".root")

		rv.writeBinToPdf(startBin, stdCmd = [folder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
		rv.plotData = True
		rv.plotTheo = False

		continue
		for b in range(startBin, stopBin):
#			if  not b == 32:
#				continue 
			rv.replaceFromRootFile("f0s.root",2)
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			intensNames = []
			argandNames = []

#			rv.writeAmplFiles(b,0,"./filesForMisha/"+sect+"_m"+str(b)+"_t"+str(tBin)+"_corrected")
#			rv.writeAmplFiles(b,1,"./filesForMisha/"+sect+"_m"+str(b)+"_t"+str(tBin)+"_uncorrect")
#			continue

			rv.writeBinToPdf(b, stdCmd = ["", folder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  folder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
	print studyList
	print cumulWeights


if __name__ == "__main__":
	main()
