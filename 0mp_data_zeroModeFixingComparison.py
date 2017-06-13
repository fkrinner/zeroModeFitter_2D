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
mPi      = 0.13957018
mK       = 0.493677

mRho     =  .77549
Grho     =  .1491

def doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rho = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 1, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"0-+0+[pi,pi]1--PiP":[rho]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData()
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
	return fitRho

def doF0phase(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]0++PiS" : (0.34, 2*mK - .04)}):
	pipiSfileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_0mp0pSigmaPiS"
	pipiSw = pc.fixedParameterization(pipiSfileName, polynomialDegree  = 0, complexPolynomial = False)
	waveModel = {"0-+0+[pi,pi]0++PiS": [pipiSw]}
	fitPiPiSshape = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitPiPiSshape.loadData()
	fitPiPiSshape.finishModelSetup()
	fitPiPiSshape.phaseFit()
	fitPiPiSshape.mode = PHASE
#	fitPiPiSshape.removeGlobalPhaseFromComa()
	fitPiPiSshape.unifyComa()
	return fitPiPiSshape

alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

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

#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_MC.root"
	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_std11.root"
	sectors          = ["0-+0+[pi,pi]0++PiS", "0-+0+[pi,pi]1--PiP"]
	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	tBins            = [tBin]
	startBin         = 12
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

	print "Starting with fixed shape f0"
	fixedShapeF0 = doFixedShapes(inFileName, sectors[:1], startBin, stopBin, tBins)
	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"

#	print "Starting with fixed shape rho"
#	fixedShapeRho = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins)
#	allMethods["fixedShapeRho"] = fixedShapeRho
#	print "Finished with fixed shape rho"

#	print "Starting with restricted rho (1 Gamma)"
#	fixedShapeRho1G = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho - Grho, mRho+Grho)})
#	allMethods["fixedShapeRho1G"] = fixedShapeRho1G
#	print "Finished with restricted rho (1 Gamma)"

#	print "Starting with restricted rho (2 Gammas)"
#	fixedShapeRho2G = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho -2*Grho, mRho+2*Grho)})
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
#	fitRho1G = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho - Grho, mRho+Grho)})
#	allMethods["fitRho1G"] = fitRho1G
#	print "Finished with fitting restricted rho (1 Gamma)"

#	print "Starting with fitting restricted rho (2 Gammas)"
#	fitRho2G = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho -2*Grho, mRho+2*Grho)})
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
	with  modernplotting.toolkit.PdfWriter("compositions_0mp_data"+str(tBin)+".pdf") as pdfOutput:
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
#	return 
	weightedSum = weightedParametersSum(evals, selfEvals, params)
	for i,m in enumerate(studyList):
		evl = sumUp(allMethods[m].evaluateZeroModeParametersForMode(weightedSum)).real
		diff = (evl - selfEvals[m])/selfEvals[m]
		evl2 = sumUp(allMethods[m].evaluateZeroModeParametersForMode(resovedWeightedSum)).real
		diff2 = (evl2 - selfEvals[m])/selfEvals[m]

		print m,diff,";:;:;;>>>??"
		hist.SetBinContent(len(studyList)+1, i+1, diff)
		hist.SetBinContent(len(studyList)+2, i+1, diff2)



	axolotl = []
	for i,study in enumerate(studyList):
		axolotl.append(shortlabels[study])
#		axolotl.append(alphabet[i])

	with modernplotting.toolkit.PdfWriter("studies_0mp_data"+str(tBin)+".pdf") as pdfOutput:
		plot = style.getPlot2D()
		plot.axes.get_xaxis().set_ticks([(i + 0.5) for i in range(len(studyList)+2)])
		plot.axes.get_yaxis().set_ticks([(i + 0.5) for i in range(len(studyList))])
		studyPlotter.makeValuePlot(plot, hist)
		
		plot.axes.set_yticklabels(axolotl)
		axolotl.append(r"$\Sigma$")
		axolotl.append(r"$\Omega$")
		plot.axes.set_xticklabels(axolotl, rotation = 90)
		plot.setZlim((0.,1.))

		pdfOutput.savefigAndClose()

	with open("studies_0mp_data"+str(tBin)+".txt",'w') as out:
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
				fileName = "./collectedMethodsData/"+stu+"_"+sect+"_0mpData_"+str(bin)
				if not (sect,bin) in fileNames:
					fileNames[sect,bin] = []
				fileNames[sect,bin].append(fileName)
				rv.writeAmplFiles(bin, fileName = fileName)

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
		allMethods['fixedShapes'].removeZeroModeFromComa()
		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = allMethods['fixedShapes'].produceResultViewer(resovedWeightedSum,s, noRun = True, plotTheory = True)
		rv.writeBinToPdf(startBin, stdCmd = ["./comparisonResultsData/" + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
		for b in range(startBin, stopBin):
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			rv.writeBinToPdf(b, stdCmd = ["", "./comparisonResultsData/" + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", intensNames,  "./comparisonResultsData/" + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", argandNames])
	print studyList
	print cumulWeights


if __name__ == "__main__":
	main()
