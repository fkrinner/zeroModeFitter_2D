from waveNameClass import waveName
from rootfabi import root_open, GetKeyNames
import numpy as np
import numpy.linalg as la
from allBinsClass import allBins
from massBinClass import massBin
from resultViewerClass import resultViewer
import sys
from utils import zeroForSectors, getZeroHistBorders, renormToBinWidth
import parameterizationClasses as pc
import scipy.optimize
from removeZeroModes import removeCertainZeroModes
from fixedparameterizationPaths import getFileNameForSector

def main():
	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_MC.root"
#	sectors = ["0-+0+[pi,pi]1--PiP","0-+0+[pi,pi]0++PiS"]
#	sectors = ["1++0+[pi,pi]1--PiS","1++0+[pi,pi]0++PiP"]
	sectors = ["2-+0+[pi,pi]0++PiD","2-+0+[pi,pi]1--PiP","2-+0+[pi,pi]1--PiF","2-+0+[pi,pi]2++PiS"]
#	sectors = ["2-+0+[pi,pi]1--PiP"]
#	sectors = ["1-+1+[pi,pi]1--PiP"]
#	sectors = ["2++1+[pi,pi]1--PiD"]
#	sectors = ["0-+0+[pi,pi]1--PiP"]
#	sectors = ["2-+0+[pi,pi]2++PiS"]
#	sectors = ["1++0+[pi,pi]1--PiS","1++0+[pi,pi]0++PiP"]
#	sectors = ["1++0+[pi,pi]1--PiS"]

	sectors = [sectors[3]]

	tBin             = 0
	startBin         = 10
	stopBin          = 50
	polynomialDegree = 1
	modelMode        = "fixedShapes"
#	modelMode        = "BW"
	useSmooth        = False
	if modelMode == "BW" or modelMode == "2BW":
		rho = pc.breitWigner()
#		rho.parameters = [.77549, .1491]
		rho.parameters = [1.27549, .100]
		rhoPrime = pc.breitWigner()
		rhoPrime.parameters = [1.570, 0.144]
		firstWaveModel = [rho]
		if modelMode == "2BW":
			firstWaveModel.append(rhoPrime)
		waveModel = {0:firstWaveModel}
	elif modelMode == "omnesReal":
		omnes = pc.omnesFunctionPolynomial(polynomialswDegree, True)
		firstWaveModel = [omnes]
		waveModel = {0:firstWaveModel}
	elif modelMode == "omnesComplex":
		omnes = pc.omnesFunctionPolynomial(polynomialDegree, False)
		firstWaveModel = [omnes]
		waveModel = {0:firstWaveModel}
	elif modelMode == "fixedShapes":
		useBF       = False
		polyDegree  = 0
		polyComplex = True
		waveModel   = {}
		for s, sector in enumerate(sectors):
			model = []
			fileNames = getFileNameForSector(sector, useBF)
			for fn in fileNames:
				param = pc.fixedParameterization(fn, polynomialDegree  = polyDegree, complexPolynomial = polyComplex)
				model.append(param)
			waveModel[s] = model
	
	

	with root_open(inFileName, "READ") as inFile:
		histNames = GetKeyNames(inFile)
		histListReal = []
		histListImag = []
		histListNorm = []
		histListIndx = []
		for sector in sectors:
			realHistName = sector+"_"+str(tBin)+"_real"
			histReal = inFile.Get(realHistName)
			if not histReal:
				raise IOError("Could not get '" + realHistName + "' from '" + inFileName + "'")
			histListReal.append(histReal)

			imagHistName = sector+"_"+str(tBin)+"_imag"
			histImag = inFile.Get(imagHistName)
			if not histImag:
				raise IOError("Could not get '" + imagHistName + "' from '" + inFileName + "'")
			histListImag.append(histImag)

			normHistName = sector+"_"+str(tBin)+"_norm"
			histNorm = inFile.Get(normHistName)
			if not histNorm:
				raise IOError("Could not get '" + normHistName + "' from '" + inFileName + "'")
			histListNorm.append(histNorm)

			indexHistName = sector+"_"+str(tBin)+"_index"
			histIndx = inFile.Get(indexHistName)
			if not histIndx:
				raise IOError("Could not get '" + indexHistName + "' from '" + inFileName + "'")
			histListIndx.append(histIndx)
		histsToFill = []
		comaHists = []
		for mBin in range(startBin, stopBin):
			comaHistName = "COMA_"+str(tBin)+"_"+str(mBin)
			comaHist = inFile.Get(comaHistName)
			if not comaHist:
				raise IOError("Could not get '" + comaHistName + "' from '" + inFileName + "'")
			comaHists.append(comaHist)
		ab = allBins(startBin, stopBin, histListReal, histListImag, histListNorm, histListIndx, comaHists)


		zeroCount = 0
		zeroHistList  = []
		eigenHistList = []
		while True:
			zeroName  = "zero"+str(zeroCount)+"_"+str(tBin)
			eigenName = "eigen"+str(zeroCount)+"_"+str(tBin)
			zeroHist = inFile.Get(zeroName)
			if not zeroHist:
				break
			zeroCount += 1
			if not zeroForSectors(sectors, zeroHist.GetTitle()):
				continue
			zeroHistList.append(zeroHist)

			eigenHist = inFile.Get(eigenName)
			if eigenHist:
				eigenHistList.append(eigenHist)
		if (not len(eigenHistList) == 0) and (not len(eigenHistList) == len(zeroHistList)):
			raise ValueError("Number of eigenvalue histograms does not match, but is also nonzero")
		removeCertainZeroModes(zeroHistList, eigenHistList)
		for zeroHist in zeroHistList:
			borders = getZeroHistBorders(zeroHist) 
			ab.addZeroMode(borders, zeroHist)
#		ab.rotateToPhaseOfBin(10)
		ab.removeZeroModeFromComa()

#		ab.unifyComa()

		if not useSmooth:
			ab.initChi2(waveModel)
	
			totalPars = []
			for k in waveModel:
				for f in waveModel[k]:
					totalPars += f.parameters
			if not len(totalPars) == 0:
				scipy.optimize.minimize(ab.chi2, totalPars)
			for k in waveModel:
				for f in waveModel[k]:
					print f.parameters

			chi2, params = ab.chi2(returnParameters = True)
			paramsZ      = ab.linearizeZeroModeParameters(params)
			ab.setTheoryFromOwnFunctions(params, True)
			print "The final chi2 =",chi2
		else:
			A,B,C   =  ab.getSmoothnessABC()
			paramsZ = -np.dot(la.inv(A + np.transpose(A)), B)
	

		intenses = []
		reals    = []
		imags    = []
		intensD  = []
		realsD   = []
		imagsD   = []
		intensT  = []
		realsT   = []
		imagsT   = []

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

			ID = rh.Clone()
			ID.Reset()
			intensD.append(ID)
	
			rD = rh.Clone()
			rD.Reset()
			realsD.append(rD)

			iD = rh.Clone()
			iD.Reset()
			imagsD.append(iD)

			IT = rh.Clone()
			IT.Reset()
			intensT.append(IT)

			rT = rh.Clone()
			rT.Reset()
			realsT.append(rT)

			iT = rh.Clone()
			iT.Reset()
			imagsT.append(iT)

		ab.fillHistograms(paramsZ, intenses)
		ab.fillHistograms(paramsZ, reals, mode = "REAL")
		ab.fillHistograms(paramsZ, imags, mode = "IMAG")
		zeroP = [0.]*len(paramsZ)
		ab.fillHistograms(zeroP, intensD)
		ab.fillHistograms(zeroP, realsD, mode = "REAL")
                ab.fillHistograms(zeroP, imagsD, mode = "IMAG")

		if not useSmooth:
			ab.fillHistograms(zeroP, intensT, mode = "INTENSTHEO")
			ab.fillHistograms(zeroP, realsT,  mode = "REALTHEO"  )
        	        ab.fillHistograms(zeroP, imagsT,  mode = "IMAGTHEO"  )

		
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

			allIsZero = True
			for binX in range(intensT[i].GetNbinsX()):
				for binY in range(intensT[i].GetNbinsY()):
					if not intensT[i].GetBinContent(binX+1, binY+1) == 0.:
						allIsZero = False
						break
				if not allIsZero:
					break # # # # # # # # # # # # # # # # # 
			if not allIsZero:
				rv = resultViewer([intenses[i], intensD[i], intensT[i]],[reals[i], realsD[i], realsT[i]],[imags[i], imagsD[i], imagsT[i]])
			else:
					rv = resultViewer([intenses[i], intensD[i]],[reals[i], realsD[i]],[imags[i], imagsD[i]])
			rv.run()

if __name__ == "__main__":
	main()
