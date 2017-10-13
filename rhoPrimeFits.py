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
from globalDefinitions import referenceWave, mPi, mK,mRho,Grho,mRhoPrime,GrhoPrime,mF0,g1,g2,mF2,GF2,Pr0
import opp_data_zeroModeFixingComparison as opp

def doFitRhoPrime(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {},referenceWave = ""):
	rhoMass  = ptc.parameter( mRho, "rhoMass" )
	rhoWidth = ptc.parameter( Grho , "rhoWidth")
	rhoPrimeMass  = ptc.parameter( mRhoPrime, "rhoPrimeMass" )
	rhoPrimeWidth = ptc.parameter( GrhoPrime, "rhoPrimeWidth")
	ratio         = ptc.parameter( .9       , "reImRatio"    )
	rho      = ptc.relativisticBreitWigner([rhoMass,rhoWidth], mPi, mPi, mPi, 1, 0, False)
	rhoPrime = ptc.relativisticBreitWignerRatio([rhoPrimeMass,rhoPrimeWidth,ratio], mPi, mPi, mPi, 1, 0, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"1++0+[pi,pi]1--PiS":[rho, rhoPrime]}, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitRho.loadData(referenceWave = referenceWave)
	fitRho.finishModelSetup()
	fitRho.fitShapeParameters()
	fitRho.calculateNonShapeParameters()
	fitRho.mode = AMPL
#	fitRho.removeGlobalPhaseFromComa()
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
	if len(sys.argv) > 2:
		study = sys.argv[2]
		studyAdder = "_"+study
	else:
		study = "std11"
		studyAdder = ""
	print "Study: "+study

	inFileName = fileNameMap[study]
	sectors    = ["1++0+[pi,pi]0++PiP", "1++0+[pi,pi]1--PiS"]
	tBins      = [tBin]
	startBin   = 30
	stopBin    = 50
#	startBin   = 25
#	stopBin    = 28


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
	fixedShapeF0 = opp.doFixedShapes(inFileName, sectors[:1], startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"

	print "Starting with fixed shapes"
	fixedShapes = opp.doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"

	print "Starting with fitting rho"
	fitRho = opp.doFitRho(inFileName, sectors, startBin, stopBin, tBins,sectorRangeMap = {"1++0+[pi,pi]0++PiP" : (0., 1.1)}, referenceWave = referenceWave)
	allMethods["fitRho"] = fitRho
	print "Finished with fitting rho"

#	print "Starting with fitting rho'"
#	fitRhoPrime = doFitRhoPrime(inFileName, sectors, startBin, stopBin, tBins, referenceWave = referenceWave)
#	allMethods["fitRhoPrime"] = fitRhoPrime
#	print "Finished with fitting rho'"

	ndfs   = {}
	params = {}
	for m in allMethods:
		ndfs[m]   =  allMethods[m].getNDFforMode()
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

	fitRhoPrime = doFitRhoPrime(inFileName, sectors[1:], startBin, stopBin, tBins, referenceWave = referenceWave)
	startPars   = [mRho,Grho,mRhoPrime,GrhoPrime,1.]
	x,err       = fitRhoPrime.fitShapeParametersForBinRange(startPars,[0],range(stopBin-startBin), zeroModeParameters = resolvedWeightedSum)
	print x
	fitRhoPrime.setShapeParameters(x,err,resolvedWeightedSum)
	s = 0 # only one sector??
	rv = fitRhoPrime.produceResultViewer(resolvedWeightedSum,s, noRun = True, plotTheory = True)
	for b in range(startBin, stopBin):
		intensName = "./rhoPrimeFits/intens_"+str(b)+"_"+str(tBin)+".pdf"
		argandName = "./rhoPrimeFits/argand_"+str(b)+"_"+str(tBin)+".pdf"
		rv.writeBinToPdf(b, stdCmd = ["", intensName, [],  argandName, []])
	return	

if __name__ == "__main__":
	main()
