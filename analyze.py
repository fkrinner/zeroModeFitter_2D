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


def main(rhoFileName = ""):
	checkLaTeX()
# # # #  Monte-Carlo results
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_MC.root"
#	inFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/DpPiPiPi_forReleaseNote.root"
#	inFileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/ppppppppp/DpPiPiPi_versionInPaper.root"
# # # # Real-data results
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_three0pp.root"
	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/std11_withIntegral.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_bigger1pp.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_exotic.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_bigger2pp.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_bigger2mp.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_3pp.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_4pp.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_4mp.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_6mp.root"
#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_exoticMC.root"

	fourPPamplFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/ampls_4++1+rhoPiG.dat"

#	inFileName = "/nfs/mds/user/fkrinner/extensiveFreedIsobarStudies/results_std11.root"
#	sectors = ["Dp[pi,pi]0++PiS","Dp[pi,pi]1--PiP"]
#	sectors = ["Dp[pi,pi]0++PiS"]
#	sectors = ["Dp[pi,pi]1--PiP"]
# # # # # Std11 Definitions
	std11_0mp0p = ["0-+0+[pi,pi]0++PiS", "0-+0+[pi,pi]1--PiP"]
	std11_1pp0p = ["1++0+[pi,pi]0++PiP","1++0+[pi,pi]1--PiS"]
	std11_1pp1p = ["1++1+[pi,pi]1--PiS"]
	std11_2pp1p = ["2++1+[pi,pi]1--PiD"]
	std11_2mp0p = ["2-+0+[pi,pi]0++PiD","2-+0+[pi,pi]1--PiP","2-+0+[pi,pi]1--PiF","2-+0+[pi,pi]2++PiS"]
	std11_2mp1p = ["2-+1+[pi,pi]1--PiP"]

# # # # # Bigger 1++ definitions
#	sectors = std11_1pp0p + ["1++0+[pi,pi]1--PiD", "1++0+[pi,pi]2++PiP"]
# # # # # Exotic (1-+)
	sector_exotic = ["1-+1+[pi,pi]1--PiP"]
# # # # # Bigger 2++ definitions
#	sectors = std11_2pp1p + ["2++1+[pi,pi]2++PiP"]
# # # # # Bigger 2-+ definitions
#	sectors = std11_2mp0p + ["2-+0+[pi,pi]2++PiD"]
#	sectors = std11_2mp1p + ["2-+1+[pi,pi]2++PiS"]
# # # # # 3++
#	sectors = ["3++0+[pi,pi]1--PiD", "3++0+[pi,pi]2++PiP", "3++0+[pi,pi]3--PiS"]
#	sectors = ["3++0+[pi,pi]2++PiP"]
# # # # # 4++ 
#	sectors = ["4++1+[pi,pi]1--PiG", "4++1+[pi,pi]2++PiF"]
#	sectors = ["4++1+[pi,pi]1--PiG"]
#	sectors = ["4++1+[pi,pi]2++PiF"]
# # # # # 4-+
#	sectors = ["4-+0+[pi,pi]1--PiF"]
# # # # # 6-+
#	sectors = ["6-+0+[pi,pi]1--PiH"]

#	sectors = std11_0mp0p
#	sectors = sector_exotic
#	sectors = std11_1pp0p
#	sectors = std11_1pp1p
#	sectors = std11_2mp0p
#	sectors = std11_2mp1p
	sectors = std11_2pp1p

	referenceWave = "4-+0+rhoPiF"

#	sectors = std11_0mp0p[:1]

	doSpecialOneBinFit = -34 # negative values turn it off

	sectorUseMap = { # Defines, if for the given sector a theory curve will be used
		"0-+0+[pi,pi]0++PiS" : True,
		"0-+0+[pi,pi]1--PiP" : True,
		"1++0+[pi,pi]0++PiP" : True, 
		"1++0+[pi,pi]1--PiS" : True, 
		"2-+0+[pi,pi]0++PiD" : True, 
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

	tBin             = 2
	startBin         = 10
	stopBin          = 50
#	startBin         = 34
#	stopBin          = 35
#	startBin         = 25
#	stopBin          = 50
	polynomialDegree = 0
	modelMode        = "fixedShapes"
#	modelMode        = "fixedRhoExotic"
#	modelMode        = "pipiS"
#	modelMode        = "none"
#	modelMode        = "simpleBW_Dp"
#	modelMode        = "BW"
#	modelMode        = "explicitRhoFile"
	useSmooth        = False
	fitShape         = True
	phaseFit         = False
	produceTotals = False
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

	elif modelMode == "fixedRhoExotic":
		pdgMass  = 0.7690
		pdgWidth = 0.1509

		massRho       = ptc.parameter(pdgMass +0.0 , "rhoMass" )
		widthRho      = ptc.parameter(pdgWidth+0.0 , "rhoWidth")
		massRho.lock  = True
		widthRho.lock = True
		rho           = ptc.breitWigner([massRho, widthRho])
		waveModel     = {"1-+1+[pi,pi]1--PiP": [rho]}

	elif modelMode == "simpleBW_Dp":
		theWarning = """
		Keep in mind the two-step procedure: 
		-	fix with the shifted
		-	c&p the parameters
		-	set the parameters by hand 
		-	run with the non-shifted "truth" parameters
		"""
		print theWarning

		f0Mass   = ptc.parameter(.98 ,"f0Mass" )
		f0Width  = ptc.parameter(.1  , "f0Width")

#		rhoMass  = ptc.parameter(.77  + .2 + .5, "rhoMass")
#		rhoWidth = ptc.parameter(.16  + .1 + .1, "rhoWidth")

		rhoMass  = ptc.parameter(.77 + .2 + .3 , "rhoMass")
		rhoWidth = ptc.parameter(.16 + .1 + .1, "rhoWidth")


		rhoPMass  = ptc.parameter(1.4, "rhoMass")
		rhoPWidth = ptc.parameter(0.2, "rhoWidth")


#		f0Mass   = ptc.parameter(1.4 , "f0Mass"  )
#		f0Width  = ptc.parameter( .1 , "f0Width" )
#		rhoMass  = ptc.parameter( .77, "rhoMass" )
#		rhoWidth = ptc.parameter( .16, "rhoWidth")

#		f0Mass   = ptc.parameter(1., "f0Mass" )
#		f0Width  = ptc.parameter(.11 ,"f0Width")
#		rhoMass  = ptc.parameter(.75, "rhoMass")
#		rhoWidparamsFarIntens.pdfth = ptc.parameter(.18, "rhoWidth")


		f0Mass.lock   = True; f0Width.lock  = True; rhoMass.lock  = True; rhoWidth.lock = True
		rhoPMass.lock = True; rhoPWidth.lock = True


		f0  = ptc.breitWigner([f0Mass, f0Width])


		rho  = ptc.breitWigner([rhoMass, rhoWidth])
		rhoP = ptc.breitWigner([rhoPMass, rhoPWidth])
#		f2 = pc.breitWigner()
#		f2.setParameters([1.27, .2 ])
#		waveModel = {"Dp[pi,pi]0++PiS" : [f0], "Dp[pi,pi]1--PiP" : [rho], "Dp[pi,pi]2++PiD" : [f2]}
		waveModel = {"Dp[pi,pi]0++PiS" : [f0], "Dp[pi,pi]1--PiP" : [rho]}

	elif modelMode == "fixedShapes":
		useBF       = True
		merge0pp    = False
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

	elif modelMode == "pipiS":
		pipiSfileName = "/nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/anything/zeroModes/bwAmplitudes_noBF/amp_0mp0pSigmaPiS"
		pipiSw = pc.fixedParameterization(pipiSfileName, polynomialDegree  = 0, complexPolynomial = False)
		waveModel = {"0-+0+[pi,pi]0++PiS": [pipiSw]}


	elif modelMode == "none":
		waveModel = {}
		for sector in sectors:
			waveModel[sector] = []

	else:
		raise RuntimeError("modelMode: '" + modelMode + "' is unknown")
	
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
		histsToFill  = []
		comaHists    = []
		intHistsReal = []
		intHistsImag = []
		for mBin in range(startBin, stopBin):
			comaHistName = "COMA_"+str(tBin)+"_"+str(mBin)
			comaHist = inFile.Get(comaHistName)

			if produceTotals:
				realHistName = "INTEGRAL_r_"+str(tBin)+"_"+str(mBin)
				realHist = inFile.Get(realHistName)
				if not realHist:
					raise IOError("Could not get '" + realHistName + "' from '" + inFileName + "'")
				intHistsReal.append(realHist)
				imagHistName = "INTEGRAL_i_"+str(tBin)+"_"+str(mBin)
				imagHist = inFile.Get(imagHistName)
				if not imagHist:
					raise IOError("Could not get '" + imagHistName + "' from '" + inFileName + "'")
				intHistsImag.append(imagHist)

			if not comaHist:
				raise IOError("Could not get '" + comaHistName + "' from '" + inFileName + "'")
			comaHists.append(comaHist)

		if not referenceWave == "":
			refHistReal  = inFile.Get(referenceWave+"_"+str(tBin)+"_real")
			if not refHistReal:
				raise IOError("Could not get '" + referenceWave+"_"+str(tBin)+"_real" + "' from '" + self.inFileName + "'")
			refHistImag  = inFile.Get(referenceWave+"_"+str(tBin)+"_imag")
			if not refHistImag:
				raise IOError("Could not get '" + referenceWave+"_"+str(tBin)+"_imag" + "' from '" + self.inFileName + "'")
			refHistIndex = inFile.Get(referenceWave+"_"+str(tBin)+"_index")
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

		rotateToFourPPampls = True
		if rotateToFourPPampls:
			fourPPampls = loadAmplsTM(fourPPamplFileName)
			ab.removePhases(fourPPampls[tBin][startBin:stopBin])

#		ab.rotateToPhaseOfBin(10)
#		ab.removeZeroModeFromComa()
#		ab.removeGlobalPhaseFromComa()
#		ab.unifyComa()

		if phaseFit:
			ab.setMassRanges(sectorRangeMap)
			for mb in ab.massBins:
				mb.setZeroTheory()
			from random import random
			ab.initChi2(waveModel)
			nBin = startBin
			leBin = ab.massBins[nBin - startBin]
			nTries = 10
			mnn = float("inf")
#			for tr in range(nTries):
#				print "at",tr,'/',nTries
#				pars  = [random() for _ in range(leBin.nParAll())]
##				print leBin.phaseChi2(pars)
#				res = scipy.optimize.minimize(leBin.phaseChi2, pars)
#				if res.fun < mnn:
#					print "Improved from",mnn,"to",res.fun,"in try",tr
#					mnn = res.fun
#					bestPar = res.x[:]
#			print bestPar, "bestPar"
			paramsPhaseChi2 = [148.57310258,  143.66171657,  116.67030827 ,  -6.8412118]
			paramsAmplChi2  = [-52.51465293, -10.3576874 ,   6.69282325,  64.28970961]
			print "----------------====----------------------------"
			leBin.phaseChi2(paramsPhaseChi2)
			print "----------------====----------------------------"
			leBin.phaseChi2(paramsAmplChi2)
			print "----------------====----------------------------"
#			paramsZ = [bestPar[0], bestPar[1]] * 1 # Very, very bad hack...
#			paramsZ = paramsPhaseChi2[:2]
			ptu = paramsAmplChi2

			leBin.phaseChi2(ptu)
			paramsZ = ptu[:2]

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

		if not useSmooth and not phaseFit and not modelMode == "none":
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

			for k in waveModel:
				for f in waveModel[k]:
					print "function parameters",f.getParameters()
			chi2, params = ab.chi2(shapePars,returnParameters = True)
			errs = ab.getNonShapeUncertainties(shapePars)
#			print "pars", len(params[0])
#			print "errs",len(errs[0])
#			print params,"params"
			print "fitted = [",
			for i in range(1, len(params[0])/2):
				if not i == 1:
					print ',',
				print params[0][2*i], '+'+ str(params[0][2*i+1])+'j',
			print ']'
			paramsZ      = ab.linearizeZeroModeParameters(params)
			ab.setTheoryFromOwnFunctions(params, True)
			print "The final chi2 =",chi2
		elif not phaseFit:
			A,B,C   =  ab.getSmoothnessABC()
			paramsZ = -np.dot(la.inv(A + np.transpose(A)), B)

		if modelMode == 'none':
			ab.initChi2(waveModel)
			params = []
			for mb in ab.massBins:
				params.append([0.]*mb.nFunc*2)
				
			ab.setTheoryFromOwnFunctions(params, True)

#		for i in range(len(paramsZ)):
#			paramsZ[i] = 0.

#		ab.writeZeroModeCoefficients(paramsZ, "zeroModeCorrections_std11", str(tBin))

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

		print paramsZ

#		print paramsZ
#		paramsZ = np.asarray([ 329.73899894 , 150.29589973]) # For the wrong parameterizations in the D decay
#		paramsZ = np.asarray([ -49.77267381,    2.93486152]) # For the P-only for in 0-+0+ at bin 25
#		ab.removeGlobalPhaseFromComa()
		ab.removeZeroModeFromComa()

		if produceTotals:
			hists = []
			for m in range(len(sectors)):
				hists.append(get3PiHistogram(sectors[m] + "_t"+ str(tBin)))
				ab.fillTotal(paramsZ, hists)

			with root_open("totals_noZeroModeWavesFromStd11.root", "UPDATE"):
				for h in hists:
					h.Write()
			return


		ab.fillHistograms(paramsZ, intenses            )
		ab.fillHistograms(paramsZ, reals,  mode = REAL )
		ab.fillHistograms(paramsZ, imags,  mode = IMAG )
		ab.fillHistograms(paramsZ, phases, mode = PHASE)
		ab.fillHistograms(paramsZ, correl, mode = REIMCORRELATION)

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
 
			noRun = False
			if not allIsZero:
				rv = resultViewer([intenses[i], intensD[i], intensT[i]],[reals[i], realsD[i], realsT[i]],[imags[i], imagsD[i], imagsT[i]], [phases[i], phasesD[i], phasesT[i]], startBin = startBin, reImCorrel = ric, noRun = noRun)
			else:
				rv = resultViewer([intenses[i], intensD[i]],[reals[i], realsD[i]],[imags[i], imagsD[i]],[phases[i], phasesD[i]], startBin = startBin, reImCorrel = ric, noRun = noRun)
			rv.titleRight = getProperWaveName(sectors[i])
#			rv.titleRight = r"$0^{-+}0^+[\pi\pi]_{1^{--}}\pi P$"
#			rv.tString    = ""
			rv.tString    = getProperDataSet(inFileName, tBin)
#			rv.tString    = "Monte Carlo"
			rv.printLiminary = False
			rv.scaleTo = "maxCorrData"
#			rv.plotCorr   = False
#			rv.plotTheo   = False
			rv.run()
#			for j in range(startBin, stopBin):
#				fileName = "./uncorrectedPlots/" + sectors[i] + "_bin" + str(j) + "_tBin"+str(tBin) + ".pdf"
#				rv.writeBinToPdf(j, stdCmd = ["",fileName,[],"",[]])
if __name__ == "__main__":
	main()
