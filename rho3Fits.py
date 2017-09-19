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

mRho3 = 1.6888
Grho3 = 0.161

def doFitRho3(inFileName,zeroFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	rho3Mass  = ptc.parameter( mRho3, "rhoMass" )
	rho3Width = ptc.parameter( Grho3, "rhoWidth")
	rho3 = ptc.relativisticBreitWigner([rho3Mass,rho3Width], mPi, mPi, mPi, 3, 0, False)
	fitRho = amplitudeAnalysis(inFileName, sectors, {"3++0+[pi,pi]3--PiS":[rho3]}, startBin, stopBin, tBins, sectorRangeMap =  {"3++0+[pi,pi]3--PiS":(1.4,2.0)}, zeroFileName = zeroFileName)
	fitRho.loadData(referenceWave = "4-+0+rhoPiF")
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

	inFileName   = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_3pp.root"
	zeroFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_3pp.root"
	sectors          = ["3++0+[pi,pi]3--PiS"]
	tBin = int(sys.argv[1])
	if tBin < 0 or tBin > 3:
		raise ValueError("Invalid t' bin: " + str(tBin))

	tBins            = [tBin]

	startBin         = 32
	stopBin          = 50

	methodBinRanges = {
	                   "fitF2"        : (22, 50),
	                   "fitRho3"      : (32, 50),
	                   "fitBothF2"    : (22, 50),
	                   "fitRhoF"      : (22, 50),
	                   "fixedShapeF2" : (22, 50)}
#	methodBinRanges = {} # Override here 


	print "Starting with fitting rho3"
	fitRho3 = doFitRho3(inFileName,zeroFileName, sectors, startBin, stopBin, tBins)
	print "Finished with fitting rho3"
	s = 0 # only one sector??
	rv = fitRho3.produceResultViewer(fitRho3.getZeroModeParametersForMode(),s, noRun = True, plotTheory = True)
	for b in range(startBin, stopBin):
		intensName = "./rho3fits/intens_"+str(b)+"_"+str(tBin)+".pdf"
		argandName = "./rho3fits/argand_"+str(b)+"_"+str(tBin)+".pdf"
		rv.writeBinToPdf(b, stdCmd = ["", intensName, [],  argandName, []])
	return 

if __name__ == "__main__":
	main()
