import matplotlib
matplotlib.use("Agg")

from analysisClass import amplitudeAnalysis
import parameterTrackingParameterizations as ptc
import parameterizationClasses as pc
from fixedparameterizationPaths import getFileNameForSector
from modes import PHASE, AMPL, SMOOTH, NONE
from utils import sumUp, weightedSum
import sys
import pyRootPwa

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

def doF0phase(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]0++PiS" : (0., 2*mK - .04)}):
	pipiSfileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_0mp0pSigmaPiS"
	pipiSw = pc.fixedParameterization(pipiSfileName, polynomialDegree  = 0, complexPolynomial = False)
	waveModel = {"0-+0+[pi,pi]0++PiS": [pipiSw]}
	fitPiPiSshape = amplitudeAnalysis(inFileName, sectors, waveModel, startBin, stopBin, tBins, sectorRangeMap = sectorRangeMap)
	fitPiPiSshape.loadData()
	fitPiPiSshape.finishModelSetup()
	fitPiPiSshape.phaseFit()
	fitPiPiSshape.mode = PHASE
	fitPiPiSshape.removeGlobalPhaseFromComa()
	return fitPiPiSshape

alphapet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def doFixedShapes(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {}):
	waveModel = {}
	for n,sector in enumerate(sectors):
		model = []
		fileNames = getFileNameForSector(sector, False, True)
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
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_MC.root"
	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_std11.root"
	sectors          = ["0-+0+[pi,pi]0++PiS", "0-+0+[pi,pi]1--PiP"]
	tBins            = [0]
	startBin         = 10
	stopBin          = 50

	allMethods       = {}

	print "Starting with fixed shape f0"
	fixedShapeF0 = doFixedShapes(inFileName, sectors[:1], startBin, stopBin, tBins)
	allMethods["fixedShapeF0"] = fixedShapeF0
	print "Finished with fixed shape f0"

	print "Starting with fixed shape rho"
	fixedShapeRho = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins)
	allMethods["fixedShapeRho"] = fixedShapeRho
	print "Finished with fixed shape rho"

	print "Starting with restricted rho (1 Gamma)"
	fixedShapeRho1G = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho - Grho, mRho+Grho)})
	allMethods["fixedShapeRho1G"] = fixedShapeRho1G
	print "Finished with restricted rho (1 Gamma)"

	print "Starting with restricted rho (2 Gammas)"
	fixedShapeRho2G = doFixedShapes(inFileName, sectors[1:], startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho -2*Grho, mRho+2*Grho)})
	allMethods["fixedShapeRho2G"] = fixedShapeRho2G
	print "Finished with restricted rho (2 Gammas)"

	print "Starting with fixed shapes"
	fixedShapes = doFixedShapes(inFileName, sectors, startBin, stopBin, tBins)
	allMethods["fixedShapes"] = fixedShapes
	print "Finished with fixed shapes"

	print "Starting with phase"
	fitPiPiSshape = doF0phase(inFileName, sectors[:1], startBin, stopBin, tBins)
	allMethods["pipiS"] = fitPiPiSshape
	print "Finished with phase"

	print "Starting with fitting restricted rho (1 Gamma)"
	fitRho1G = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho - Grho, mRho+Grho)})
	allMethods["fitRho1G"] = fitRho1G
	print "Finished with fitting restricted rho (1 Gamma)"

	print "Starting with fitting restricted rho (2 Gammas)"
	fitRho2G = doFitRho(inFileName, sectors, startBin, stopBin, tBins, sectorRangeMap = {"0-+0+[pi,pi]1--PiP":(mRho -2*Grho, mRho+2*Grho)})
	allMethods["fitRho2G"] = fitRho2G
	print "Finished with fitting restricted rho (2 Gammas)"

	if stopBin - startBin > 1:
		print "Starting with smooth"
		smooth = doSmooth(inFileName, sectors, startBin, stopBin, tBins)
		allMethods["smooth"] = smooth
		print "Finished with smooth"

	fileNames = {}

	for stu in allMethods:
		print "Writing for '" + stu + "'"
		for s, sect in enumerate(allMethods[stu].sectors):
			rv = allMethods[stu].produceResultViewer(allMethods[stu].getZeroModeParametersForMode(),s, noRun = True)
			for bin in range(startBin, stopBin):
				fileName = "./collectedMethods/"+stu+"_"+sect+"_0mpMC_"+str(bin)
				if not (sect,bin) in fileNames:
					fileNames[sect,bin] = []
				fileNames[sect,bin].append(fileName)
				rv.writeAmplFiles(bin, fileName = fileName)
	studyList = []
	for m in allMethods:
		studyList.append(m)
	studyList.sort()
	params    = {}
	selfEvals = {}
	for m in allMethods:
		params[m] = allMethods[m].getZeroModeParametersForMode()
		selfEvals[m] = sumUp(allMethods[m].evaluateZeroModeParametersForMode(params[m])).real

	hist = pyRootPwa.ROOT.TH2D("hist","hist", len(params)+1, 0, len(params)+1, len(params), 0, len(params))

	evals = {}
	for i,m in enumerate(studyList):
		print "-------------------------------"
		for j,n in enumerate(studyList):
			evl = sumUp(allMethods[n].evaluateZeroModeParametersForMode(params[m])).real
			evals[n,m] = evl
			diff = (evl-selfEvals[n])/selfEvals[n]
			print m,'in',n,":",diff
			hist.SetBinContent(i+1, j+1, diff)

	weightedSum = weightedParametersSum(evals, selfEvals, params)
	for i,m in enumerate(studyList):
		evl = sumUp(allMethods[m].evaluateZeroModeParametersForMode(weightedSum)).real
		diff = (evl - selfEvals[m])/selfEvals[m]
		print m,diff,";:;:;;>>>??"
		hist.SetBinContent(len(studyList)+1, i+1, diff)



	style = modernplotting.mpplot.PlotterStyle()
#	style.p2dColorMap = 'ocean_r'
#	style.p2dColorMap = 'YlOrRd'
	style.p2dColorMap = 'Reds'
	axolotl = []
	for i in range(len(studyList)):
		axolotl.append(alphapet[i])
#		axolotl.append('')

	with modernplotting.toolkit.PdfWriter("studies_0mp.pdf") as pdfOutput:
		plot = style.getPlot2D()
		plot.axes.get_xaxis().set_ticks([])
		plot.axes.get_yaxis().set_ticks([])
		studyPlotter.makeValuePlot(plot, hist)

		mpsp.plotCustomAxesTicks(plot, axolotl,direction='y', offset = .5)
		axolotl.append(r"$\Sigma$")
		mpsp.plotCustomAxesTicks(plot, axolotl, offset = .5)

		plot.axes.get_xaxis().set_ticks([])
		plot.axes.get_yaxis().set_ticks([])
		plot.setZlim((0.,1.))

		pdfOutput.savefigAndClose()

	for s, sect in enumerate(allMethods['fixedShapes'].sectors):
		allMethods['fixedShapes'].removeGlobalPhaseFromComa()
		rv = allMethods['fixedShapes'].produceResultViewer(weightedSum,s, noRun = True)
		rv.writeBinToPdf(startBin, stdCmd = ["./comparisonResults/" + sect + "_MC_2D.pdf", "", [], "", []])
		for b in range(startBin, stopBin):
			intensNames = [name+".intens" for name in fileNames[sect,b]]
			argandNames = [name+".argand" for name in fileNames[sect,b]]
			rv.writeBinToPdf(b, stdCmd = ["", "./comparisonResults/" + sect + "_MC_intens_"+str(b)+".pdf", intensNames,  "./comparisonResults/" + sect + "_MC_argand_"+str(b)+".pdf", argandNames])	
	print studyList


if __name__ == "__main__":
	main()
