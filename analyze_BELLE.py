#!/usr/bin/python
from waveNameClass import waveName
from rootfabi import root_open, GetKeyNames
import numpy as np
import numpy.linalg as la
from allBinsClass import allBins
from massBinClass import massBin
from modes import REAL, IMAG, PHASE, INTENS, INTENSNORM, INTENSTHEO, REALTHEO, IMAGTHEO, PHASETHEO, REIMCORRELATION
from resultViewerClass import resultViewer
import os, sys
from utils import zeroForSectors, getZeroHistBorders, renormToBinWidth, checkLaTeX, changeReferenceWave
import parameterizationClasses as pc
import parameterTrackingParameterizations as ptc
import scipy.optimize
from removeZeroModes import removeCertainZeroModes
from fixedparameterizationPaths import getFileNameForSector
from LaTeX_strings import getProperWaveName, getProperDataSet
from utils import get3PiHistogram, loadAmplsTM

sys.path.append("/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp")
from analyzeBELLE_version2 import parseCmdLine
from getBranchFileEnding   import getBranchFileEnding

def main(rhoFileName = ""):
	checkLaTeX()

	freeMap, freeString = parseCmdLine(sys.argv)

	inFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/DdecayResults_"+freeString+"_"+getBranchFileEnding()+".root"

	zeroFileName = inFileName

	sectors = []
	allSectors = ["KpiSright","KpiPright", "KpiDright", "KpiSwrong","KpiPwrong","KpiDwrong", "piPiS","piPiP","piPiD"]

	fixMap = [False, False, False, False, False, False, False, False, False]

	for i in range(9):
		if freeString[i] == '1':
			if not fixMap[i]:
				sectors.append(allSectors[i])

#	sectors = sectors[:1]

	sectorRangeMap   = {}

	tBin             = 0
	startBin         = 34
	stopBin          = 35

	modelMode        = "simpleBW_Dp"
	if modelMode == "simpleBW_Dp":

		Kmass  = ptc.parameter(1.414, "KMass")
		Kwidth = ptc.parameter(.232, "KWidth")
		Kmass.lock = True; Kwidth.lock = True
		Kstar = ptc.breitWigner([Kmass, Kwidth])

		K892Mass  = ptc.parameter(.89166, "K892Mass")
		K892Width = ptc.parameter(.0508, "K902Width")
		K892Mass.lock = True; K892Width.lock = True
		K892 = ptc.breitWigner([K892Mass, K892Width])

		f0Mass   = ptc.parameter(.98 ,"f0Mass" )
		f0Width  = ptc.parameter(.1  , "f0Width")
		f0Mass.lock   = True; f0Width.lock  = True; 
		f0    = ptc.breitWigner([f0Mass, f0Width])

		rhoMass  = ptc.parameter(.77526, "rhoMass")
		rhoWidth = ptc.parameter(.1478, "rhoWidth")
		rhoMass.lock  = True; rhoWidth.lock = True
		rho  = ptc.breitWigner([rhoMass, rhoWidth])

		f2Mass   = ptc.parameter( 1.2751,"f2Mass" )
		f2Width  = ptc.parameter( .1851 , "f2Width")
		f2Mass.lock   = True; f2Width.lock  = True; 
		f2    = ptc.breitWigner([f2Mass, f2Width])
		
		waveModel = {"KpiSright":[Kstar], "KpiPright":[K892],"KpiSwrong":[Kstar],"KpiPwrong":[K892],"piPiS":[f0],"piPiP":[rho],"piPiD":[f2]}

	uniCOMA = False
	with root_open(inFileName, "READ") as inFile:
		histNames = GetKeyNames(inFile)
		histListReal = []
		histListImag = []
		histListNorm = []
		histListIndx = []
		nAmplMax     = 0
		for sector in sectors:
			realHistName = sector+"_"+str(tBin)+"_real"
			histReal = inFile.Get(realHistName)
			if not histReal:
				raise IOError("Could not get '" + realHistName + "' from '" + inFileName + "'")
			histReal.SetDirectory(0)
			histListReal.append(histReal)

			imagHistName = sector+"_"+str(tBin)+"_imag"
			histImag = inFile.Get(imagHistName)
			if not histImag:
				raise IOError("Could not get '" + imagHistName + "' from '" + inFileName + "'")
			histImag.SetDirectory(0)
			histListImag.append(histImag)

			normHistName = sector+"_"+str(tBin)+"_norm"
			histNorm = inFile.Get(normHistName)
			if not histNorm:
				raise IOError("Could not get '" + normHistName + "' from '" + inFileName + "'")
			histNorm.SetDirectory(0)
			histListNorm.append(histNorm)

			indexHistName = sector+"_"+str(tBin)+"_index"
			histIndx = inFile.Get(indexHistName)
			if not histIndx:
				raise IOError("Could not get '" + indexHistName + "' from '" + inFileName + "'")
			histIndx.SetDirectory(0)
			nAmplMax = max(nAmplMax, int(2*histIndx.GetMaximum()))
			histListIndx.append(histIndx)
		histsToFill  = []
		comaHists    = []
		intHistsReal = []
		intHistsImag = []
		for mBin in range(startBin, stopBin):
			comaHistName = "COMA_"+str(tBin)+"_"+str(mBin)
			comaHist = inFile.Get(comaHistName)
			if not comaHist:
				raise IOError("Could not load '" + comaHistName + "'")

			comaHist.SetDirectory(0)
			comaHists.append(comaHist)

		ab = allBins(startBin, stopBin, histListReal, histListImag, histListNorm, histListIndx, comaHists, intHistsReal, intHistsImag)

	with root_open(zeroFileName, "READ") as inFile:
		zeroCount = 0
		zeroHistList  = []
		eigenHistList = []
		while True:
			zeroName  = "zero"+str(zeroCount)+"_"+str(tBin)
			eigenName = "eigen"+str(zeroCount)+"_"+str(tBin)
			zeroHist  = inFile.Get(zeroName)
			if not zeroHist:
				break
			zeroHist.SetDirectory(0)
			print "Adding zero-mode" 
			zeroCount += 1
			if not zeroForSectors(sectors, zeroHist.GetTitle()):
				continue
			zeroHistList.append(zeroHist)

			eigenHist = inFile.Get(eigenName)
			eigenHist.SetDirectory(0)
			if eigenHist:
				eigenHistList.append(eigenHist)
		if (not len(eigenHistList) == 0) and (not len(eigenHistList) == len(zeroHistList)):
			raise ValueError("Number of eigenvalue histograms does not match, but is also nonzero")
		removeCertainZeroModes(zeroHistList, eigenHistList)
		for zeroHist in zeroHistList:
			borders = getZeroHistBorders(zeroHist) 
			ab.addZeroMode(borders, zeroHist)

		rotateToFourPPampls = False
		if rotateToFourPPampls:
			fourPPampls = loadAmplsTM(fourPPamplFileName)
			ab.removePhases(fourPPampls[tBin][startBin:stopBin])

#		ab.rotateToPhaseOfBin(10)
		zeroModeComaVal = 100. #float(sys.argv[1])
#		zeroModeComaVal = 100000.
#		zeroModeComaVal = "unchanged"
#		ab.removeZeroModeFromComa()
#		ab.addComaValueForZeroMode(zeroModeComaVal)
#		ab.removeGlobalPhaseFromComa()
		ab.unifyComa()

		ab.initChi2(waveModel)
		ab.setMassRanges(sectorRangeMap)
		totalPars = []
		for k in waveModel:
			for f in waveModel[k]:
				totalPars += f.getParameters()
		shapePars = []
		if not len(totalPars) == 0 and fitShape:
			res = scipy.optimize.minimize(ab.chi2, totalPars)
			hi = res.hess_inv
			print "m0 = ",res.x[0],"+-",(2*hi[0,0])**.5
			print "G0 = ",res.x[1],"+-",(2*hi[1,1])**.5
			print "Giving a Chi2 of:",res.fun
			shapePars = res.x

		chi2, params = ab.chi2(shapePars,returnParameters = True)
		errs         = ab.getNonShapeUncertainties(shapePars)
		paramsZ      = ab.linearizeZeroModeParameters(params)
		zmPar        = ab.getNonShapeParameters()
		ab.setTheoryFromOwnFunctions(params, True)

		intenses = []
		reals    = []
		imags    = []
		correl   = []
		phases   = []
		intensD  = []
		realsD   = []
		imagsD   = []
		phasesD  = []
		intensT  = []
		realsT   = []
		imagsT   = []
		phasesT  = []

		for rh in histListReal:
			Ih = rh.Clone()
			Ih.Reset()
			intenses.append(Ih)

			realH = rh.Clone()
			realH.Reset()
			reals.append(realH)
	
			imagH = rh.Clone()
			imagH.Reset()
			imags.append(imagH)

			reImCorrH = rh.Clone()
			reImCorrH.Reset()
			correl.append(reImCorrH)

			phaseH = rh.Clone()
			phaseH.Reset()
			phases.append(phaseH)

			ID = rh.Clone()
			ID.Reset()
			intensD.append(ID)
	
			rD = rh.Clone()
			rD.Reset()
			realsD.append(rD)

			iD = rh.Clone()
			iD.Reset()
			imagsD.append(iD)

			pD = rh.Clone()
			pD.Reset()
			phasesD.append(pD)

			IT = rh.Clone()
			IT.Reset()
			intensT.append(IT)

			rT = rh.Clone()
			rT.Reset()
			realsT.append(rT)

			iT = rh.Clone()
			iT.Reset()
			imagsT.append(iT)

			pT = rh.Clone()
			pT.Reset()
			phasesT.append(pT)

		zeroP = [0.]*len(paramsZ)
#		paramsZ = zeroP

		ab.fillHistograms(paramsZ, intenses            )
		ab.fillHistograms(paramsZ, reals,  mode = REAL )
		ab.fillHistograms(paramsZ, imags,  mode = IMAG )
		ab.fillHistograms(paramsZ, phases, mode = PHASE)
		ab.fillHistograms(paramsZ, correl, mode = REIMCORRELATION)

		ab.fillHistograms(zeroP, intensD              )
		ab.fillHistograms(zeroP, realsD,  mode = REAL )
                ab.fillHistograms(zeroP, imagsD,  mode = IMAG )
                ab.fillHistograms(zeroP, phasesD, mode = PHASE)

		ab.fillHistograms(zeroP, intensT,  mode = INTENSTHEO)
		ab.fillHistograms(zeroP, realsT,   mode = REALTHEO  )
       	        ab.fillHistograms(zeroP, imagsT,   mode = IMAGTHEO  )
       	        ab.fillHistograms(zeroP, phasesT,  mode = PHASETHEO )

		for i in range(len(histListReal)):
			renormToBinWidth(intenses[i]  )
			renormToBinWidth(intensD[i]   )
			renormToBinWidth(intensT[i]   )
			renormToBinWidth(reals[i] , .5)
			renormToBinWidth(realsD[i], .5)
			renormToBinWidth(realsT[i], .5)
			renormToBinWidth(imags[i] , .5)
			renormToBinWidth(imagsD[i], .5)
			renormToBinWidth(imagsT[i], .5)
			renormToBinWidth(correl[i],   )
			allIsZero = True
			for binX in range(intensT[i].GetNbinsX()):
				for binY in range(intensT[i].GetNbinsY()):
					if not intensT[i].GetBinContent(binX+1, binY+1) == 0.:
						allIsZero = False
						break
				if not allIsZero:
					break # # # # # # # # # # # # # # # # #
			ric = correl[i]
 
			noRun = False
			if not allIsZero:
				rv = resultViewer([intenses[i], intensD[i], intensT[i]],[reals[i], realsD[i], realsT[i]],[imags[i], imagsD[i], imagsT[i]], [phases[i], phasesD[i], phasesT[i]], startBin = startBin, reImCorrel = ric, noRun = noRun)
			else:
				rv = resultViewer([intenses[i], intensD[i]],[reals[i], realsD[i]],[imags[i], imagsD[i]],[phases[i], phasesD[i]], startBin = startBin, reImCorrel = ric, noRun = noRun)
			rv.titleRight    = getProperWaveName(sectors[i])
			rv.tString       = getProperDataSet(inFileName, tBin)
	                rv.titleFontSize = 11
			rv.overrideMassString = ""

			rv.printLiminary = False
			rv.topMarginIntens = 1.4
			rv.run()
#			for j in range(startBin, stopBin):
#				fileName = "./uncorrectedPlots/" + sectors[i] + "_bin" + str(j) + "_tBin"+str(tBin) + ".pdf"
#				rv.writeBinToPdf(j, stdCmd = ["",fileName,[],"",[]])
if __name__ == "__main__":
	main()
