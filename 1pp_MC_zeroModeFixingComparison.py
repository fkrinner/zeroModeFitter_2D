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

import consistencyUtils as cu

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

from estimateErrors import estimateErrors2

import scipy.optimize
mPi      = 0.13957018
mK       = 0.493677

#mRho     =  .77549
#Grho     =  .1491

mRho = .7
Grho = 0.1

def doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	rhoMass  = ptc.parameter( mRho , "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	binning  = [0.278, 0.32, 0.36, 0.4 , 0.44, 0.48, 0.52, 0.56, 0.6 , 0.64, 0.68, 0.7 , 0.72,
	            0.74 , 0.76, 0.78, 0.8 , 0.82, 0.84, 0.86, 0.88, 0.9 , 0.92, 0.96, 1.0 , 1.04, 
	            1.08 , 1.12, 1.16, 1.2 , 1.24, 1.28, 1.32, 1.36, 1.4 , 1.44, 1.48, 1.52, 1.56, 
	            1.6  , 1.64, 1.68, 1.72, 1.76, 1.8 , 1.84, 1.88, 1.92, 1.96, 2.0 , 2.04, 2.08, 
	            2.12 , 2.16, 2.2 , 2.24, 2.28]

#	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 0, False)
	rho = ptc.integratedRelativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 0, binning, intensWeight = True, fitPr = False, reweightInverseBW = True)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"1++0+[pi,pi]1--PiS":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData()
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doF0phase(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"1++0+[pi,pi]0++PiP" : (0.34, 2*mK - .04)}):
	pipiSfileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_1pp0pSigmaPiP"
	pipiSw = pc.fixedParameterization(pipiSfileName, polynomialDegree  = 0, complexPolynomial = False)
	waveModel = {"1++0+[pi,pi]0++PiP": [pipiSw]}
	fitPiPiSshape = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitPiPiSshape.loadData()
	fitPiPiSshape.finishModelSetup()
	fitPiPiSshape.phaseFit()
	fitPiPiSshape.mode = PHASE
#	fitPiPiSshape.removeGlobalPhaseFromComa()
	fitPiPiSshape.unifyComa()
	return fitPiPiSshape

def doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
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
	fixedShapes.loadData()
	fixedShapes.finishModelSetup()
	fixedShapes.fitShapeParameters()
	fixedShapes.calculateNonShapeParameters()
	fixedShapes.mode = AMPL
#	fixedShapes.removeGlobalPhaseFromComa()
	return fixedShapes

def doSmooth(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	waveModel = {}
	smooth = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
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

	style.titleRight = r"$1^{++}0^+$"
	style.titleLeft  = r"Monte-Carlo"


	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_MC.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_std11.root"
	sectors          = ["1++0+[pi,pi]0++PiP", "1++0+[pi,pi]1--PiS"]
	tBins            = [0]
#	startBin         = 19
#	stopBin          = 50
	startBin = 12
	stopBin  = 50

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


	print "Starting with fixed shape f0"
	fixedShapeF0 = doFixedShapes(inFileName, sectors[:1], startBin, stopBin, tBins)
	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"

#	print "Starting with fixed shape rho"
#	fixedShapeRho = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins)
#	allMethods["fixedShapeRho"] = fixedShapeRho
#	print "Finished with fixed shape rho"

#	print "Starting with restricted rho (1 Gamma)"
#	fixedShapeRho1G = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, sectorRangeMap = {"1++0+[pi,pi]1--PiS":(mRho - Grho, mRho+Grho)})
#	allMethods["fixedShapeRho1G"] = fixedShapeRho1G
#	print "Finished with restricted rho (1 Gamma)"

#	print "Starting with restricted rho (2 Gammas)"
#	fixedShapeRho2G = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, sectorRangeMap = {"1++0+[pi,pi]1--PiS":(mRho -2*Grho, mRho+2*Grho)})
#	allMethods["fixedShapeRho2G"] = fixedShapeRho2G
#	print "Finished with restricted rho (2 Gammas)"

	print "Starting with fixed shapes"
	fixedShapes = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"

#	print "Starting with phase"
#	fitPiPiSshape = doF0phase(inFileName, sectors[:1], startBin, stopBin, tBins)
#	allMethods["pipiS"] = fitPiPiSshape
#	print "Finished with phase"

	print "Starting with fitting rho"
	fitRho = doFitRho(inFileName, sectors, startBin, stopBin, tBins)
	allMethods["fitRho"] = fitRho
	print "Finished with fitting rho"

#	print "Starting with fitting restricted rho (1 Gamma)"
#	fitRho1G = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"1++0+[pi,pi]1--PiS":(mRho - Grho, mRho+Grho)})
#	allMethods["fitRho1G"] = fitRho1G
#	print "Finished with fitting restricted rho (1 Gamma)"

#	print "Starting with fitting restricted rho (2 Gammas)"
#	fitRho2G = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"1++0+[pi,pi]1--PiS":(mRho -2*Grho, mRho+2*Grho)})
#	allMethods["fitRho2G"] = fitRho2G
#	print "Finished with fitting restricted rho (2 Gammas)"

	if stopBin - startBin > 1:
		print "Starting with smooth"
		smooth = doSmooth(inFileName, sectors, startBin, stopBin, tBins)
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
	with  modernplotting.toolkit.PdfWriter("compositions_1pp.pdf") as pdfOutput:
		plot = style.getPlot1D()
		for m in comps:
			line  = [0.]*len(comps[m][0])
			xAxis = [ .5 + 0.04*(startBin + i) for i in range(len(comps[m][0]))]
			break
		count = 0
		for m in comps:
			print m

			newLine = line[:]
			for i in range(len(comps[m][0])):
				newLine[i] += comps[m][0][i]
			plot.axes.fill_between(xAxis, line, newLine, facecolor = modernplotting.colors.makeColorLighter(modernplotting.colors.colorScheme.blue, 0.1*count))
			count += 1
			line = newLine
		plot.setYlim(0.,1.)
		plot.setXlim(xAxis[0], xAxis[-1])
		pdfOutput.savefigAndClose()
	studyList = []
	for m in allMethods:
		studyList.append(m)
	studyList.sort()

	selfEvals = {}
	for m in allMethods:
		selfEvals[m] = sumUp(allMethods[m].evaluateResolvedZeroModeParametersForMode(params[m])).real

	hist = pyRootPwa.ROOT.TH2D("hist","hist", len(params)+2, 0, len(params)+2, len(params), 0, len(params))

	cumulWeights = {}
	resovedWeightedSum = [[]] # Assumes one t' bin
	for i in range(stopBin - startBin):
		dim = len(params[m][0][i])
		prrs = [0.] * dim
		for m in params:
			weight = comps[m][0][i]
			if not m in cumulWeights:
				cumulWeights[m] = 0.
			cumulWeights[m] += weight
			for j in range(dim):
#				print dim,i,j,params[m][0]
#				print prrs
				prrs[j] += weight * params[m][0][i][j]
		resovedWeightedSum[0].append(prrs)

########################### THE BLOCK #########################################
	if False:
		print "Begin the second fit"
		fitRho.setExternalZMP(resovedWeightedSum)
		startVal = fitRho.fitParameters[:]
		startVal[0] += 0.2
		startVal[1] += 0.2
		res = scipy.optimize.minimize(fitRho.callAmplChi2WithFixedZeroMode,startVal)
		hi        = res.hess_inv
	
		parsaren  = []
		errs      = []

		for i in range(len(res.x)):
			errs.append((2.*hi[i,i])**.5)
			parsaren.append(res.x[i])
		print "errs before", errs
		errs = estimateErrors2(fitRho.callAmplChi2WithFixedZeroMode, parsaren, errs)

		pars = res.x[:]
		for i in range(len(res.x)):
			print "-/-/-/-/-/-/-/-/-/-",i,"-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-"
			pars[i] += errs[i]
			print fitRho.callAmplChi2WithFixedZeroMode(pars) - res.fun
			pars[i] -= 2*errs[i]
			print fitRho.callAmplChi2WithFixedZeroMode(pars) - res.fun
			pars[i] += errs[i]
		print "-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-"

		for i in range(2):
			print res.x[i],"+-", errs[i]

		if res.x[0] < 1.05:
			print "Compare with rho parameters"
			print "sigM =", (abs(res.x[0]-0.7690)/errs[0])
			print "sigG =", (abs(res.x[1]-0.1509)/errs[1])
		else:
			print "Compare with f2 parameters"
			print "sigM =", (abs(res.x[0]-1.2751)/errs[0])
			print "sigG =", (abs(res.x[1]-0.1851)/errs[1])
		return 
################################################################################
	evals = {}
	for i,m in enumerate(studyList):
#		print "-------------------------------"
		for j,n in enumerate(studyList):
			evl = sumUp(allMethods[n].evaluateResolvedZeroModeParametersForMode(params[m])).real
			evals[n,m] = evl
			diff = (evl-selfEvals[n])/selfEvals[n]
			
#			allMethods["fixedShapes"].removeZeroModeFromComa()
#			print "------------------------------------IN---------------------------------"
#			print params[m], params[n]
#			diff = sumUp(allMethods["fixedShapes"].compareTwoZeroModeCorrections(params[m], params[n]))
#			print diff
#			print "------------------------------------OUT---------------------------------"
#			print m,'in',n,":",diff
			hist.SetBinContent(i+1, j+1, diff)
	weightedSum = weightedParametersSum(evals, selfEvals, params)
	for i,m in enumerate(studyList):
		evl = sumUp(allMethods[m].evaluateZeroModeParametersForMode(cloneZeros(weightedSum))).real
		diff = (evl - selfEvals[m])/selfEvals[m]
		evl2 = sumUp(allMethods[m].evaluateZeroModeParametersForMode(resovedWeightedSum)).real
		diff2 = (evl2 - selfEvals[m])/selfEvals[m]

		print m,diff,";:;:;;>>>??"
		hist.SetBinContent(len(studyList)+1, i+1, diff)
		hist.SetBinContent(len(studyList)+2, i+1, diff2)

	axolotl = []
	for i,study in enumerate(studyList):
		axolotl.append(shortlabels[study])

	style.titleRight = r"$1^{++}0^+$"
	style.titleLeft  = r"Monte Carlo"

	with modernplotting.toolkit.PdfWriter("studies_1pp.pdf") as pdfOutput:
		plot = style.getPlot2D()
		plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(len(studyList)+2)])
		plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(len(studyList))])
		studyPlotter.makeValuePlot(plot, hist)
		
		plot.axes.set_yticklabels(axolotl)
#		axolotl.append(r"$\Sigma$")
		axolotl.append(r"$\vec 0$")
		axolotl.append(r"$\Omega$")
		plot.axes.set_xticklabels(axolotl, rotation = 90)
		plot.setZlim((0.,1.))

		pdfOutput.savefigAndClose()

	with open("studies_1pp.txt",'w') as out:
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
				fileName = "./collectedMethods/"+stu+"_"+sect+"_1ppMC_"+str(bin)
				if not (sect,bin) in fileNames:
					fileNames[sect,bin] = []
				fileNames[sect,bin].append(fileName)
				rv.writeAmplFiles(bin, fileName = fileName)

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
		allMethods['fixedShapes'].removeZeroModeFromComa()
		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
#		print "Achtung::: cloneZeros() is aktiv"
		rv = allMethods['fixedShapes'].produceResultViewer(resovedWeightedSum,s, noRun = True, plotTheory = True)
		rv.writeBinToPdf(startBin, stdCmd = ["./comparisonResults/" + sect + "_MC_2D.pdf", "", [], "", []])

#		continue # remove me

		for b in range(startBin, stopBin):
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			rv.writeBinToPdf(b, stdCmd = ["", "./comparisonResults/" + sect + "_MC_intens_"+str(b)+".pdf", intensNames,  "./comparisonResults/" + sect + "_MC_argand_"+str(b)+".pdf", argandNames])
	print studyList
	print cumulWeights


if __name__ == "__main__":
	main()
