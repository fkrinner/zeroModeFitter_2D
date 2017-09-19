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

import LaTeX_strings
import consistencyUtils as cu
from studyFileNames import fileNameMap

import modernplotting.mpplot
import modernplotting.toolkit
import modernplotting.specialPlots as mpsp
import studyPlotter

import opp_data_zeroModeFixingComparison as opp

mPi      = 0.13957018
mK       = 0.493677

mRho     =  .77549
Grho     =  .1491

#mRho = .7
#Grho = 0.1
Pr0  = 0.1973

mRhoPrime = 1.465
GrhoPrime =  .400


def main():
	checkLaTeX()
	style = modernplotting.mpplot.PlotterStyle()
#	style.p2dColorMap = 'ocean_r'
#	style.p2dColorMap = 'YlOrRd'
	style.p2dColorMap = 'Reds'

	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))
	if len(sys.argv) > 2:
		study = sys.argv[2]
		studyAdder = "_"+study
	else:
		study = "std11"
		studyAdder = ""
	print "Study: "+study

	inFileName = fileNameMap[study]
	sectors          = ["1++0+[pi,pi]0++PiP", "1++0+[pi,pi]1--PiS"]
	tBins            = [tBin]
	startBin = 13
	stopBin  = 50
#	startBin = 25
#	stopBin  = 28


	allMethods       = {}
	methodStrings    = {}
	shortlabels      = {  "fixedShapeF0"    : r"$\text{fix}_{f_0}^{~}$",
	                      "fixedShapeRho"   : r"$\text{fix}_\rho^{~}$",
	                      "fixedShapeRho1G" : r"$\text{fix}_\rho^{1\Gamma}$",
	                      "fixedShapeRho2G" : r"$\text{fix}_\rho^{2\Gamma}$",
	                      "fixedShapes"     : r"$\text{fix}_\text{all}^{~}$",
	                      "pipiS"           : r"$\phi_{[\pi\pi]_S}^{~}$",
	                      "fitRho"          : r"$\text{fit}_\rho^{~}$",
	                      "fitRhoPrime"     : r"$\text{fit}_{\rho^\prime}^{~}$",
	                      "fitRho1G"        : r"$\text{fit}_\rho^{1\Gamma}$",
	                      "fitRho2G"        : r"$\text{fit}_\rho^{2\Gamma}$",
	                      "smooth"          : r"smooth"}

	print "Starting with fixed shape f0"
	fixedShapeF0 = opp.doFixedShapes(inFileName, sectors[:1], startBin, stopBin, tBins)
	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"

	print "Starting with fixed shapes"
	fixedShapes = opp.doFixedShapes(inFileName, sectors, startBin, stopBin, tBins)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"

	print "Starting with fitting rho"
	fitRho = opp.doFitRho(inFileName, sectors, startBin, stopBin, tBins,sectorRangeMap = {"1++0+[pi,pi]0++PiP" : (0., 1.1)})
	allMethods["fitRho"] = fitRho
	print "Finished with fitting rho"

	print "Starting with fitting rho'"
	fitRhoPrime = opp.doFitRhoPrime(inFileName, sectors, startBin, stopBin, tBins)
	allMethods["fitRhoPrime"] = fitRhoPrime
	print "Finished with fitting rho'"

	ndfs   = {}
	params = {}
	for m in allMethods:
		ndfs[m]=  allMethods[m].getNDFforMode()
		params[m] = allMethods[m].getZeroModeParametersForMode()
		print m,sumUp(allMethods[m].evaluateZeroModeParametersForMode(params[m])).real/ndfs[m]
	diffs = cu.getmBinResolvedDiffs(allMethods)
	comps = cu.getCompositions(diffs)
	studyList = []
	for m in allMethods:
		studyList.append(m)
	studyList.sort()

	selfEvals = {}
	for m in allMethods:
		selfEvals[m] = sumUp(allMethods[m].evaluateResolvedZeroModeParametersForMode(params[m])).real

	hist = pyRootPwa.ROOT.TH2D("hist","hist", len(params)+2, 0, len(params)+2, len(params), 0, len(params))

	cumulWeights = {}
	resolvedWeightedSum = [[]] # Assumes one t' bin
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
		resolvedWeightedSum[0].append(prrs)

	evals = {}
	for i,m in enumerate(studyList):
#		print "-------------------------------"
		for j,n in enumerate(studyList):
			evl = sumUp(allMethods[n].evaluateResolvedZeroModeParametersForMode(params[m])).real
			evals[n,m] = evl
			diff = (evl-selfEvals[n])/selfEvals[n]
			
			hist.SetBinContent(i+1, j+1, diff)
	weightedSum = opp.weightedParametersSum(evals, selfEvals, params)
	for i,m in enumerate(studyList):
		evl = sumUp(allMethods[m].evaluateZeroModeParametersForMode(cloneZeros(weightedSum))).real
		diff = (evl - selfEvals[m])/selfEvals[m]
		evl2 = sumUp(allMethods[m].evaluateZeroModeParametersForMode(resolvedWeightedSum)).real
		diff2 = (evl2 - selfEvals[m])/selfEvals[m]

		print m,diff,";:;:;;>>>??"
		hist.SetBinContent(len(studyList)+1, i+1, diff)
		hist.SetBinContent(len(studyList)+2, i+1, diff2)


	axolotl = []
	for i,study in enumerate(studyList):
		axolotl.append(shortlabels[study])
#		axolotl.append(alphabet[i])

	style.titleRight = r"$1^{++}0^+$"
	style.titleLeft  = LaTeX_strings.tBins[tBin]


	with modernplotting.toolkit.PdfWriter("studies_1pp_data"+str(tBin)+studyAdder+".pdf") as pdfOutput:
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

	with open("studies_1pp_data"+str(tBin)+studyAdder+".txt",'w') as out:
		for axl in axolotl:
			out.write(axl + ' ')
		out.write("\n")
		for i in range(hist.GetNbinsX()):
			for j in range(hist.GetNbinsY()):
				out.write(str(hist.GetBinContent(i+1, j+1)) + ' ')
			out.write('\n')

	folder = "./rhoPrimeFits/"

	for s, sect in enumerate(fitRhoPrime.sectors):
		if "0++" in sect:
			continue
		fitRhoPrime.removeZeroModeFromComa()
		fitRhoPrime.removeGlobalPhaseFromComa()
		rv = fitRhoPrime.produceResultViewer(resolvedWeightedSum,s, noRun = True, plotTheory = True)
		rv.writeBinToPdf(startBin, stdCmd = [folder + sect + "_data_2D_"+str(tBin)+".pdf", "", [], "", []])
		for b in range(startBin, stopBin):
			rv.writeBinToPdf(b, stdCmd = ["", folder + sect + "_data_intens_"+str(b)+"_"+str(tBin)+".pdf", [],  folder + sect + "_data_argand_"+str(b)+"_"+str(tBin)+".pdf", []])
	return	

if __name__ == "__main__":
	main()
