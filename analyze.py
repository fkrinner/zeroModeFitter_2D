from waveNameClass import waveName
from rootfabi import root_open, GetKeyNames
import numpy as np
import numpy.linalg as la
from allBinsClass import allBins
from massBinClass import massBin
from modes import REAL, IMAG, PHASE, INTENS, INTENSNORM, INTENSTHEO, REALTHEO, IMAGTHEO, PHASETHEO, REIMCORRELATION
from resultViewerClass import resultViewer
import os, sys
from utils import zeroForSectors, getZeroHistBorders, renormToBinWidth
import parameterizationClasses as pc
import scipy.optimize
from removeZeroModes import removeCertainZeroModes
from fixedparameterizationPaths import getFileNameForSector


def main(rhoFileName = ""):
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_MC.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_3pp.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_bigger2pp.root"
	inFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/DpPiPiPi.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_std11.root"
	sectors = ["Dp[pi,pi]0++PiS","Dp[pi,pi]1--PiP"]
#	sectors = ["Dp[pi,pi]0++PiS","Dp[pi,pi]2++PiD"]
#	sectors = ["Dp[pi,pi]1--PiP","Dp[pi,pi]2++PiD"]
#	sectors = ["Dp[pi,pi]0++PiS","Dp[pi,pi]1--PiP","Dp[pi,pi]2++PiD"]
#	sectors = ["0-+0+[pi,pi]0++PiS", "0-+0+[pi,pi]1--PiP"]
#	sectors = ["1++0+[pi,pi]0++PiP","1++0+[pi,pi]1--PiS"]
#	sectors = ["2-+0+[pi,pi]0++PiD","2-+0+[pi,pi]1--PiP","2-+0+[pi,pi]1--PiF","2-+0+[pi,pi]2++PiS"]
#	sectors = ["2-+0+[pi,pi]1--PiP"]
#	sectors = ["2-+0+[pi,pi]0++PiD"]
#	sectors = ["1-+1+[pi,pi]1--PiP"]
#	sectors = ["1++0+[pi,pi]0++PiP"]
#	sectors = ["2++1+[pi,pi]1--PiD"]
#	sectors = ["0-+0+[pi,pi]0++PiS"]
#	sectors = ["0-+0+[pi,pi]1--PiP"]
#	sectors = ["0-+0+[pi,pi]0++PiS", "0-+0+[pi,pi]1--PiP"]
#	sectors = ["2-+0+[pi,pi]2++PiS"]
#	sectors = ["2-+1+[pi,pi]1--PiP"]
#	sectors = ["1++0+[pi,pi]1--PiS","1++0+[pi,pi]0++PiP"]
#	sectors = ["1++1+[pi,pi]1--PiS"]
#	sectors = ["1++0+[pi,pi]1--PiS"]
#	sectors = ["3++0+[pi,pi]3--PiS"]
#	sectors = ["3++0+[pi,pi]2++PiP"]
#	sectors = ["3++0+[pi,pi]1--PiD"]
#	sectors = ["3++0+[pi,pi]3--PiS","3++0+[pi,pi]2++PiP","3++0+[pi,pi]1--PiD"]
#	sectors = ["4++1+[pi,pi]1--PiG", "4++1+[pi,pi]2++PiF"]
#	sectors = ["2+q+1+[pi,pi]1--PiD","2++1+[pi,pi]2++PiP"]

#	sectors = [sectors[1]]

	doSpecialOneBinFit = -15 # negative values turn it off

	sectorUseMap = { # Defines, if for the given sector a theory curve will be used
		"0-+0+[pi,pi]0++PiS" : False,
		"0-+0+[pi,pi]1--PiP" : True,
		"1++0+[pi,pi]0++PiP" : False, 
		"1++0+[pi,pi]1--PiS" : True, 
		"2-+0+[pi,pi]0++PiD" : False, 
		"2-+0+[pi,pi]1--PiP" : True, 
		"2-+0+[pi,pi]1--PiF" : True, 
		"2-+0+[pi,pi]2++PiS" : True
	}

	sectorRangeMap = {}
	for sector in sectors:
		if '1--' in sector:
			pass
#			sectorRangeMap[sector] = (0.55,1.00)
		if '0++' in sector:
			pass
#			sectorRangeMap[sector] = (0.,.94)

	tBin             = 0
	startBin         = 34
	stopBin          = 35
#	startBin         = 15
#	stopBin          = 50
	polynomialDegree = 0
#	modelMode        = "
	modelMode        = "simpleBW_Dp"
#	modelMode        = "BW"
#	modelMode        = "explicitRhoFile"
	useSmooth        = False
	fitShape         = True
	phaseFit         = False

	if modelMode == "BW" or modelMode == "2BW":
#		rho = pc.breitWigner()
		rho = pc.rpwaBreitWignerInt(0.13957018,0.13957018,0.13957018,1,0)
		rhoParameters = [.77549, .1491]
#		rhoParameters = [1.2, .2]
		rho.setParameters(rhoParameters)
#		rho.weightIntens = True
#		rhoParameters = [1.27549, .100]
		rhoPrime = pc.breitWigner()
		rhoPrimeParameters = [1.570, 0.144]
		rhoPrime.setParameters(rhoPrimeParameters)
		rhoModel = [rho]
		if modelMode == "2BW":
			rhoModel.append(rhoPrime)
		waveModel = {}
		for sector in sectors:
			if '1--' in sector:
				waveModel[sector] = rhoModel
			if "3--" in sector:
				waveModel[sector] = rhoModel

	elif modelMode == "simpleBW_Dp":
		f0  = pc.breitWigner()
		f0.setParameters([.98, .1])
		rho = pc.breitWigner()
		rho.setParameters([ .77, .16])
#		f2 = pc.breitWigner()
#		f2.setParameters([1.27, .2 ])
#		waveModel = {"Dp[pi,pi]0++PiS" : [f0], "Dp[pi,pi]1--PiP" : [rho], "Dp[pi,pi]2++PiD" : [f2]}
		waveModel = {"Dp[pi,pi]0++PiS" : [f0], "Dp[pi,pi]1--PiP" : [rho]}

	elif modelMode == "fixedShapes":
		useBF       = False
		merge0pp    = True
		polyDegree  = 0
		polyComplex = True
		waveModel   = {}
		for sector in sectors:
#			if "3--" in sector:
#				rho = pc.rpwaBreitWignerInt(0.13957018,0.13957018,0.13957018,1,1)
#				rho.setParameters([1.2, .2])
#				model = [rho]
#				waveModel[sector] = model
#				continue
			if sector in sectorUseMap:
				if not sectorUseMap[sector]:
					continue
			model = []
			fileNames = getFileNameForSector(sector, useBF, merge0pp)
			print fileNames
			for fn in fileNames:
				param = pc.fixedParameterization(fn, polynomialDegree  = polyDegree, complexPolynomial = polyComplex)
				model.append(param)
			waveModel[sector] = model
#		for sector in sectors: # Ovveride with free rho parameterization
#			if '1--' in sector:
#				rho = pc.rpwaBreitWignerInt(0.13957018,0.13957018,0.13957018,1,0) # Override here
#				rho.setParameters([.77549, .1491])
#				waveModel[sector] = [rho]

	elif modelMode == "explicitRhoFile":
		if not os.path.isfile(rhoFileName):
			raise IOError("Rho file does not exist")
		param   = pc.fixedParameterization(rhoFileName, polynomialDegree  = 0, complexPolynomial = False)
		waveModel = {}
		for sector in sectors:
			if '1--' in sector:
				waveModel[sector] = [param]

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
			print "Adding zero-mode" 
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
		ab.removeGlobalPhaseFromComa()
#		ab.unifyComa()

		if phaseFit:
			for mb in ab.massBins:
				mb.setZeroTheory()
			from random import random
			ab.initChi2(waveModel)
			nBin = 25
			leBin = ab.massBins[nBin - startBin]
			pars  = [random() for _ in range(leBin.nParAll())]
			print leBin.phaseChi2(pars)
			res = scipy.optimize.minimize(leBin.phaseChi2, pars)
			print "phaseFit gives for bin",nBin 
			print "Giving a Chi2 of:",res.fun
			paramsZ = [res.x[0], res.x[1]] * 100 # Very, very bad hack...
	
		if doSpecialOneBinFit >= 0:
			specialBin = ab.massBins[doSpecialOneBinFit]
			specialBin.initChi2(waveModel)
			print "Doing special fit for: ",+ specialBin.bin3pi
			pars = [.77549, .1491]
			if modelMode == "explicitRhoFile":
				return specialBin.chi2([])
			res = scipy.optimize.minimize(specialBin.chi2, pars)
			hi = res.hess_inv
			print "m0 = ",res.x[0],"+-",hi[0,0]**.5
			print "G0 = ",res.x[1],"+-",hi[1,1]**.5
			print "Giving a Chi2 of:",res.fun

			sys.exit(0)

		if not useSmooth and not phaseFit:
			ab.initChi2(waveModel)
			ab.setMassRanges(sectorRangeMap)
			totalPars = []
			for k in waveModel:
				for f in waveModel[k]:
					totalPars += f.getParameters()
			if not len(totalPars) == 0 and fitShape:
				res = scipy.optimize.minimize(ab.chi2, totalPars)
				hi = res.hess_inv
				print "m0 = ",res.x[0],"+-",hi[0,0]**.5
				print "G0 = ",res.x[1],"+-",hi[1,1]**.5
				print "Giving a Chi2 of:",res.fun
			for k in waveModel:
				for f in waveModel[k]:
					print "function parameters",f.getParameters()
			chi2, params = ab.chi2(returnParameters = True)
			paramsZ      = ab.linearizeZeroModeParameters(params)
			ab.setTheoryFromOwnFunctions(params, True)
			print "The final chi2 =",chi2
		elif not phaseFit:
			A,B,C   =  ab.getSmoothnessABC()
			paramsZ = -np.dot(la.inv(A + np.transpose(A)), B)

#		for i in range(len(paramsZ)):
#			paramsZ[i] = 0.

		ab.writeZeroModeCoefficients(paramsZ, "zeroModeCorrections_std11", str(tBin))

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

		ab.fillHistograms(paramsZ, intenses            )
		ab.fillHistograms(paramsZ, reals,  mode = REAL )
		ab.fillHistograms(paramsZ, imags,  mode = IMAG )
		ab.fillHistograms(paramsZ, phases, mode = PHASE)
		ab.fillHistograms(paramsZ, correl, mode = REIMCORRELATION)
		zeroP = [0.]*len(paramsZ)
		ab.fillHistograms(zeroP, intensD              )
		ab.fillHistograms(zeroP, realsD,  mode = REAL )
                ab.fillHistograms(zeroP, imagsD,  mode = IMAG )
                ab.fillHistograms(zeroP, phasesD, mode = PHASE)

		if not useSmooth:
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
 
			if not allIsZero:
				rv = resultViewer([intenses[i], intensD[i], intensT[i]],[reals[i], realsD[i], realsT[i]],[imags[i], imagsD[i], imagsT[i]], [phases[i], phasesD[i], phasesT[i]], startBin = startBin, reImCorrel = ric)
			else:
				rv = resultViewer([intenses[i], intensD[i]],[reals[i], realsD[i]],[imags[i], imagsD[i]],[phases[i], phasesD[i]], startBin = startBin, reImCorrel = ric)
			rv.run()

if __name__ == "__main__":
	main()
