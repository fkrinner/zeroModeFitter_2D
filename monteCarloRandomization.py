import pyRootPwa
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
from comaUtils import writeComa

def main(rhoFileName = ""):
	checkLaTeX()
	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_exotic.root"
# # # # # Exotic (1-+)
	sectors = ["1-+1+[pi,pi]1--PiP"]

	referenceWave = "4-+0+rhoPiF"
#	referenceWave = ""

	mBin        = 32
	tBin        = 3
	startBin    = mBin
	stopBin     = mBin+1

	useBF       = False
	merge0pp    = False
	polyDegree  = 0
	polyComplex = True
	waveModel   = {}

	model = []
	fileNames = getFileNameForSector(sectors[0], useBF, merge0pp)
	for fn in fileNames:
		param = pc.fixedParameterization(fn, polynomialDegree  = polyDegree, complexPolynomial = polyComplex)
		model.append(param)
	waveModel[sectors[0]] = model

	with root_open(inFileName, "READ") as inFile:
		histNames = GetKeyNames(inFile)
		histListReal = []
		histListImag = []
		histListNorm = []
		histListIndx = []
		nAmplMax     = 0
		for sector in sectors:
			realHistName = sector+"_"+str(tBin)+"_real"
			histReal    = inFile.Get(realHistName)
			if not histReal:
				raise IOError("Could not get '" + realHistName + "' from '" + inFileName + "'")
			histReal.SetDirectory(0)
			histListReal.append(histReal)

			imagHistName = sector+"_"+str(tBin)+"_imag"
			histImag     = inFile.Get(imagHistName)
			if not histImag:
				raise IOError("Could not get '" + imagHistName + "' from '" + inFileName + "'")
			histImag.SetDirectory(0)
			histListImag.append(histImag)

			normHistName = sector+"_"+str(tBin)+"_norm"
			histNorm     = inFile.Get(normHistName)
			if not histNorm:
				raise IOError("Could not get '" + normHistName + "' from '" + inFileName + "'")
			histNorm.SetDirectory(0)
			histListNorm.append(histNorm)

			indexHistName = sector+"_"+str(tBin)+"_index"
			histIndx      = inFile.Get(indexHistName)
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
				print "Could not get '" + comaHistName + "' from '" + inFileName + "'", "USING IDENTITY"
				comaHist = pyRootPwa.ROOT.TH2D(comaHistName, comaHistName, nAmplMax, 0.,1.,nAmplMax, 0.,1.)
				for i in range(nAmplMax):
					comaHist.SetBinContent(i+1, i+1, 1.)

			comaHist.SetDirectory(0)
			comaHists.append(comaHist)

		if not referenceWave == "":
			refHistReal  = inFile.Get(referenceWave+"_"+str(tBin)+"_real")
			refHistReal.SetDirectory(0)
			if not refHistReal:
				raise IOError("Could not get '" + referenceWave+"_"+str(tBin)+"_real" + "' from '" + self.inFileName + "'")
			refHistImag  = inFile.Get(referenceWave+"_"+str(tBin)+"_imag")
			refHistImag.SetDirectory(0)
			if not refHistImag:
				raise IOError("Could not get '" + referenceWave+"_"+str(tBin)+"_imag" + "' from '" + self.inFileName + "'")
			refHistIndex = inFile.Get(referenceWave+"_"+str(tBin)+"_index")
			refHistIndex.SetDirectory(0)
			if not refHistIndex:
				raise IOError("Could not get '" + referenceWave+"_"+str(tBin)+"_index" + "' from '" + self.inFileName + "'")
			changeReferenceWave(histListReal, histListImag, histListIndx, comaHists, refHistReal, refHistImag, refHistIndex, startBin, stopBin)

		ab = allBins(startBin, stopBin, histListReal, histListImag, histListNorm, histListIndx, comaHists, intHistsReal, intHistsImag)

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

		writeSamples = False
		nRand        = 100

		rh    = histListReal[0]
		nPts  = rh.GetNbinsY()
		valsC = np.zeros((2*nPts,nRand))
		valsD = np.zeros((2*nPts,nRand))

		for r in range(nRand+1):

			ab.initChi2(waveModel)

			shapePars    = []
			chi2, params = ab.chi2(shapePars,returnParameters = True)
			errs         = ab.getNonShapeUncertainties(shapePars)
			paramsZ      = ab.linearizeZeroModeParameters(params)

			zeroP = [0.]*len(paramsZ)
			argandNameBase = "./samplingArgands/argand_<mode>_"+str(r)+".dat"

			IC = rh.Clone()
			IC.Reset()
			ab.fillHistograms(paramsZ, [IC],  mode = INTENS )
			renormToBinWidth(IC , 1.)

			rC = rh.Clone()
			rC.Reset()
			ab.fillHistograms(paramsZ, [rC],  mode = REAL )
			renormToBinWidth(rC , .5)

			iC = rh.Clone()
			iC.Reset()
			ab.fillHistograms(paramsZ, [iC],  mode = IMAG )
			renormToBinWidth(iC , .5)

			if not r == 0: # ugly, but r = 0 is unrandomized
				for i in range(nPts):
					valsC[2*i  ,r-1] = rC.GetBinContent(mBin+1,i+1)
					valsC[2*i+1,r-1] = iC.GetBinContent(mBin+1,i+1)
	
			if writeSamples:
				with open(argandNameBase.replace("<mode>", "C"), 'w') as outFile:
					for i in range(nPts):
						m = rC.GetYaxis().GetBinCenter(i)
						outFile.write(str(m))
						for h in [IC, rC, iC]:
							v = h.GetBinContent(mBin+1, i+1)
							e = h.GetBinError(mBin+1, i+1)
							if v == 0. and e == 0.:
								continue
							outFile.write(' ' + str(v) + ' ' + str(e))
						outFile.write('\n')

			ID = rh.Clone()
			ID.Reset()
			ab.fillHistograms(zeroP,   [ID],  mode = INTENS )
			renormToBinWidth(ID , 1.)

			rD = rh.Clone()
			rD.Reset()
			ab.fillHistograms(zeroP,   [rD],  mode = REAL )
			renormToBinWidth(rD , .5)

			iD = rh.Clone()
			iD.Reset()
		        ab.fillHistograms(zeroP,   [iD],  mode = IMAG )
			renormToBinWidth(iD , .5)

			if not r == 0: # ugly, but r = 0 is unrandomized
				for i in range(nPts):
					valsD[2*i  ,r-1] = rD.GetBinContent(mBin+1,i+1)
					valsD[2*i+1,r-1] = iD.GetBinContent(mBin+1,i+1)

			if writeSamples:
				with open(argandNameBase.replace("<mode>", "D"), 'w') as outFile:
					for i in range(nPts):
						m = rD.GetYaxis().GetBinCenter(i)
						outFile.write(str(m))
						for h in [ID, rD, iD]:
							v = h.GetBinContent(mBin+1, i+1)
							e = h.GetBinError(mBin+1, i+1)
							if v == 0. and e == 0.:
								continue
							outFile.write(' ' + str(v) + ' ' + str(e))
						outFile.write('\n')

			ab.randomize() # Put the randomization at the end, to have the fist run with unrandomized points

		correl = rh.Clone()
		correl.Reset()
		ab.fillHistograms(paramsZ, [correl], mode = REIMCORRELATION)
		renormToBinWidth(correl , 1.)
		if writeSamples:
			with open("correlations_D.dat",'w') as outFile:
				for i in range(correl.GetNbinsY()):
					outFile.write(str(rD.GetBinError(mBin+1, i+1)**2)+ ' ' + str(correl.GetBinContent(mBin+1, i+1)) + ' ' + str(iD.GetBinError(mBin+1, i+1)**2)+'\n')
		ab.removeZeroModeFromComa()
		correl.Reset()
		ab.fillHistograms(paramsZ, [correl], mode = REIMCORRELATION)
		rC.Reset()
		iC.Reset()
		ab.fillHistograms(paramsZ, [rC], mode = REAL)
		ab.fillHistograms(paramsZ, [iC], mode = IMAG)
		renormToBinWidth(correl , 1. )
		renormToBinWidth(rC     ,  .5)
		renormToBinWidth(iC     ,  .5)
		with open("correlations_C.dat",'w') as outFile:
			for i in range(correl.GetNbinsY()):
				outFile.write(str(rC.GetBinError(mBin+1, i+1)**2)+' '+str(correl.GetBinContent(mBin+1, i+1))+' '+str(iC.GetBinError(mBin+1,i+1)**2)+'\n')

#	normComas = True
#	writeComa(valsC, "comaC.dat", normed = normComas)
#	writeComa(valsD, "comaD.dat", normed = normComas)
	


if __name__ == "__main__":
	main()
